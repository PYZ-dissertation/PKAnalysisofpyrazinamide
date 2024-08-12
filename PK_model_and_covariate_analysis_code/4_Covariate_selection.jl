################
#Julia code for the selection of covariates based off the best fitted base model (one compartment, 2 transit)
#Adds weight in then analyses HIV as a covariate after which, sex, age, altres, cratres, and totbres are trialled in forward selection
#Time is then tried as a covariate
#Backward selection is performed on all added covariates and on IIV#The best fit model is then analysed using vpc, gof and shrinkage
# An exposure analysis is then done comparing UAC between dosing weight bands, HIV status and time point
##################

#Load packages
using Pkg, XLSX, Pumas, PumasUtilities, Random, CSV, Chain,
CategoricalArrays, Dates, DataSets, StatsBase, DataFrames, DataFramesMeta, CairoMakie, XLSX,
AlgebraOfGraphics, Serialization, Pumas.Latexify, QuantileRegressions, OptimalDesign, Plots, DataStructures,
Statistics, StatsPlots, GLM, Makie, Polynomials, HypothesisTests, StatsBase, Loess, StatsModels
using GLM: lm, @formula
using Bioequivalence.GLM: lm, @formula
using Pumas: fit
######

Pkg.status()
#Setting working directory
pwd()
cd("")

#Loading cleaned pyrazinamide data 
PYZ_data = DataFrame(XLSX.readtable("", "", infer_eltypes = true))

###################
#Initial code to get data to correct format

# Function to convert time string to Time object
function convert_to_time(time_str)
    if ismissing(time_str)
        return missing
    else
        return Time(time_str)
    end
end

# Data cleaning and type conversion
for col in names(PYZ_data)
    if col in ["TOTBRES", "CREATRES", "WT"]
        PYZ_data[!, col] = convert(Vector{Float64}, PYZ_data[!, col])
    elseif col == "DATE"
        PYZ_data[!, col] = Date.(string.(PYZ_data[!, col]))
    elseif col == "TIMEC"
        PYZ_data[!, col] = convert_to_time.(PYZ_data[!, col])
    end
end

#Creating a log concentration column
PYZ_data[!,:LOG_DV] = log.(PYZ_data.DV)

# Handling observations at the time of dosing (EVID = 1 or AMT = 1)
@rtransform! PYZ_data :DV = :EVID == 1 ? missing : :DV
@rtransform! PYZ_data :LOG_DV = :EVID == 1 ? missing : :LOG_DV
@rtransform! PYZ_data :DV = :AMT == 1 ? missing : :DV
@rtransform! PYZ_data :LOG_DV = :AMT == 1 ? missing : :LOG_DV

# Add :DATE and :TIMEC (converted to Time object) to create :CTIM
PYZ_data = @transform PYZ_data begin
    :CTIM = :DATE .+ :TIMEC
end

# Function to calculate time after dose
function reset_time_based_on_conditions(df)
    time = similar(df.CTIM, Float64)
    current_time_offset = DateTime(0)  # Initialize with a dummy date
    
    for i in 1:nrow(df)
        if df.EVID[i] == 1
            current_time_offset = df.CTIM[i]
            time[i] = 0.0
        else
            time_diff = df.CTIM[i] - current_time_offset
            time[i] = Dates.value(time_diff) / (1000 * 60 * 60)  # Convert milliseconds to hours
        end
    end
    return time
end

# Apply the function to each group and update the Time after dose (TAD) column
grouped_df = groupby(PYZ_data, :ID)
for sdf in grouped_df
    sdf.TAD = reset_time_based_on_conditions(sdf)
    PYZ_data[findall(x -> x in sdf.ID, PYZ_data.ID), :TAD] .= sdf.TAD
end

filter!((row) -> row.TAD ≤ 35, PYZ_data)#ID 27 has some values misspecified, removing those

# Create first_D_met DataFrame with first CTIM for each ID
first_D_met = combine(groupby(PYZ_data, :ID),
               nrow => :count,
               :CTIM => first => :first_date,
               )

# Join first_D_met with PYZ_data               
PYZ_data=innerjoin(PYZ_data,
               first_D_met,
               on = :ID)

# Calculate TIME from CTIM and first_date               
PYZ_data = @transform PYZ_data begin
    :TIME = (:CTIM.-:first_date)
end

#One missing time, replacing with last known time of that ID
function fill_missing_times!(data::DataFrame, time_column::Symbol)
    last_time = nothing  # Variable to store the last non-missing time
    for i in 1:nrow(data)
        if ismissing(data[i, time_column])
            data[i, time_column] = last_time
        else
            last_time = data[i, time_column]
        end
    end
end
# Apply the function to your DataFrame
fill_missing_times!(PYZ_data, :TIME)

# Convert TIME to hours
PYZ_data = @transform PYZ_data begin
    :TIME = Dates.value.(:TIME)
end
PYZ_data = @transform PYZ_data begin
    :TIME = :TIME/(1000*60*60)
end

#Sort data by ID's the time
PYZ_data = sort(PYZ_data, [:ID, :TIME])

# Function to adjust duplicate times so they are monotonic increasing
function adjust_duplicate_times(times)
    time_counts = Dict{Float64, Int}()
    adjusted_times = similar(times)
    for i in 1:length(times)
        time = times[i]
        if haskey(time_counts, time)
            time_counts[time] += 1
            adjusted_times[i] = time + (time_counts[time] - 1) * 0.00000001
        else
            time_counts[time] = 1
            adjusted_times[i] = time
        end
    end
    return adjusted_times
end

# Apply the function to each group and update the TAD column
PYZ_data = combine(groupby(PYZ_data, :ID)) do sdf
    sdf.TAD = adjust_duplicate_times(sdf.TAD)
    sdf
end

# Function to add Time point (TP) column
function add_tp_column!(df::DataFrame)
    # Create the new TP column with zeros
    df.TP = zeros(Int, nrow(df))
    
    # Get the unique IDs
    unique_ids = unique(df.ID)
    
    # Loop over each unique ID
    for id in unique_ids
        # Get the indices for this ID
        idx = findall(df.ID .== id)
        
        # Find the first occurrence of "W6" in the Time.point column for this ID
        first_w6_index = findfirst(i -> occursin("W6", df."Time.point"[i]), idx)
        
        # If there is an occurrence of "W6", update the TP column for this ID
        if first_w6_index !== nothing
            df.TP[first_w6_index:end] .= 1
        end
    end
    
    return df
end

# Function to add Time point 2 (Time to differentiate between d0, w6 and months) column
function add_tp2_column!(df::DataFrame)
    # Create the new TP2 column with the default value of 2
    df.TP2 = fill(2, nrow(df))
    
    # Update the TP2 column based on the presence of "D0" or "W6" in the Time.point values
    df.TP2[occursin.("D0", df."Time.point")] .= 0
    df.TP2[occursin.("W6", df."Time.point")] .= 1
    
    return df
end

# Add TP column to PYZ_data
PYZ_data = add_tp_column!(PYZ_data)

# Add TP2 column to PYZ_data
PYZ_data = add_tp2_column!(PYZ_data)

# Set ROUTE column to "ev" for all rows
PYZ_data.ROUTE .= "ev"

# Filter data for D0 and W6 time points
PYZ_dataD0 = filter(row -> occursin("D0", row."Time.point"), PYZ_data)
PYZ_dataW6 = filter(row -> occursin("W6", row."Time.point"), PYZ_data)

# Sort data by ID and TAD
PYZ_dataD0 = sort(PYZ_dataD0, [:ID, :TAD])
PYZ_dataW6 = sort(PYZ_dataW6, [:ID, :TAD])

#Data now in correct format, can begin covariate analysis
#####################

#Read in a Pumas model for the full data
PYZ_pk = read_pumas(
    PYZ_data;
    id = :ID,
    time = :TIME,
    amt = :AMT,
    observations = [:LOG_DV],
    evid = :EVID,
    cmt = :CMT,
    covariates= [:WT, :HIV, :SEX, :ALTRES, :CREATRES, :AGE, :TOTBRES, :TP, :TP2]
)
#Read in a Pumas model for the D0 data
PYZ_pkD0 = read_pumas(
    PYZ_dataD0;
    id = :ID,
    time = :TAD,
    amt = :AMT,
    observations = [:LOG_DV],
    evid = :EVID,
    cmt = :CMT,
    covariates= [:WT, :HIV, :SEX, :ALTRES, :CREATRES, :AGE, :TOTBRES, :TP, :TP2]
)
#Read in a Pumas model for the W6 data
PYZ_pkW6 = read_pumas(
    PYZ_dataW6;
    id = :ID,
    time = :TAD,
    amt = :AMT,
    observations = [:LOG_DV],
    evid = :EVID,
    cmt = :CMT,
    covariates= [:WT, :HIV, :SEX, :ALTRES, :CREATRES, :AGE, :TOTBRES,:TP, :TP2]
)

#########################################
#ADDING WT TO OUR BEST FIT BASE MODEL####
PYZ_WT = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =4)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init =25)
        """
        Mean transit time
        """
        tMTT ∈ RealDomain(; lower = 0, init = 1)
        """
          - ΩCL
          - ΩVc
          - ΩMTT
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [ 0.1, 0.1, 0.1, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.38)
    end
    @covariates begin
        """
        Weight (kg)
        """
        WT
    end
    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        CL = tvcl * (WT/median(PYZ_data.WT))^0.75 * exp(η[1]) 
        Vc = tvv * (WT/median(PYZ_data.WT)) * exp(η[2])
        MTT=tMTT * exp(η[3])

        ktr = 3/MTT 
        k20= CL/Vc
    end

    @dosecontrol begin
        bioav = ( ; Abs=1 * exp(η[4]))
    end

    @dynamics begin
        Abs' = -ktr*Abs
        B1' = ktr *Abs- ktr*B1
        B2' = ktr *B1- ktr*B2
        Central' = ktr*B2 - k20 * Central
    end
    @derived begin
        cp := @. log(0.0001+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p) 
    end
end
#Initial parms
pkparam_WT = init_params(PYZ_WT)
#Fit using FOCE
pkfit_WT = fit(PYZ_WT, PYZ_pk, pkparam_WT, FOCE())
#Inference
infer(pkfit_WT)
#Calculate BIC
bic(pkfit_WT)

#########
#Testing HIV as a covariate on CL
PYZ_HIV_CL = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =1)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init = 45)
        """
        Mean transit time
        """
        tMTT ∈ RealDomain(; lower = 0, init = 2)
        """
        HIV CL factor
        """
        HIV_CL ∈ RealDomain(;  init = 0.5)
        """
          - ΩCL
          - ΩVc
          - ΩMTT
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [ 0.1, 0.1, 0.1, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.38)
    end
    @covariates begin
        """
        Weight (kg)
        """
        WT
        """
        HIV
        """
        HIV
    end
    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        CL = tvcl * (WT/median(PYZ_data.WT))^0.75 * (1+ HIV_CL*HIV) * exp(η[1]) 
        Vc = tvv * (WT/median(PYZ_data.WT)) * exp(η[2])
        MTT = tMTT * exp(η[3])

        ktr= 3/MTT
        k20= CL/Vc
    end

    @dosecontrol begin
        bioav = ( ; Abs=1 * exp(η[4]))
    end

    @dynamics begin
        Abs' = -ktr*Abs
        B1' = ktr *Abs- ktr*B1
        B2' = ktr *B1- ktr*B2
        Central' = ktr*B2 - k20 * Central
    end
    @derived begin
        cp := @. log(0.00001+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p) 
    end
end
#Initial parms
pkparam_CL = init_params(PYZ_HIV_CL)
#Fit using FOCE
pkfit_HIV_CL = fit(PYZ_HIV_CL, PYZ_pk, pkparam_CL, FOCE())
#Inference
infer(pkfit_HIV_CL)
#Calculate BIC
bic(pkfit_HIV_CL)

#Testing HIV as a covariate on V
PYZ_HIV_V = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =3)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init = 35)
        """
        Mean transit time
        """
        tMTT ∈ RealDomain(; lower = 0, init = 0.5)
        """
        HIV V factor
        """
        HIV_V ∈ RealDomain(; init = 0.5)
        """
          - ΩCL
          - ΩVc
          - ΩMTT
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [ 0.1, 0.1, 0.1, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.38)
    end
    @covariates begin
        """
        Weight (kg)
        """
        WT
        """
        HIV
        """
        HIV
    end
    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        CL = tvcl * (WT/median(PYZ_data.WT))^0.75 * exp(η[1]) 
        Vc = tvv * (WT/median(PYZ_data.WT)) * (1+ HIV_V*HIV) * exp(η[2])
        MTT = tMTT * exp(η[3])
        
        ktr = 3/MTT
        k20= CL/Vc
    end

    @dosecontrol begin
        bioav = ( ; Abs=1 * exp(η[4]))
    end

    @dynamics begin
        Abs' = -ktr*Abs
        B1' = ktr *Abs- ktr*B1
        B2' = ktr *B1- ktr*B2
        Central' = ktr*B2 - k20 * Central
    end
    @derived begin
        cp := @. log(0.00001+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p) 
    end
end
#Initial parms
pkparam_V = init_params(PYZ_HIV_V )
#Fit using FOCE
pkfit_HIV_V = fit(PYZ_HIV_V , PYZ_pk, pkparam_V, FOCE())
#Inference
infer(pkfit_HIV_V)
#Calculate BIC
bic(pkfit_HIV_V)

#Testing HIV as a covariate on MTT
PYZ_HIV_MTT = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =1)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init = 45)
        """
        Mean transit time
        """
        tMTT ∈ RealDomain(; lower = 0, init = 2)
        """
        HIV MTT factor
        """
        HIV_MTT ∈ RealDomain(; lower = 0, init = 0.5)
        """
          - ΩCL
          - ΩVc
          - ΩMTT
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [ 0.1, 0.1, 0.1, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; init = 0.38)
    end
    @covariates begin
        """
        Weight (kg)
        """
        WT
        """
        HIV
        """
        HIV
    end
    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        CL = tvcl * (WT/median(PYZ_data.WT))^0.75 * exp(η[1]) 
        Vc = tvv * (WT/median(PYZ_data.WT)) *  exp(η[2])
        MTT =tMTT * (1+ HIV_MTT*HIV) * exp(η[3])

        ktr = 3/MTT
        k20= CL/Vc
    end

    @dosecontrol begin
        bioav = ( ; Abs=1 * exp(η[4]))
    end
 
    @dynamics begin
        Abs' = -ktr*Abs
        B1' = ktr *Abs- ktr*B1
        B2' = ktr *B1- ktr*B2
        Central' = ktr*B2 - k20 * Central
    end
    @derived begin
        cp := @. log(0.00001+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p) 
    end
end
#Initial parms
pkparam_MTT = init_params(PYZ_HIV_MTT )
#Fit using FOCE
pkfit_HIV_MTT = fit(PYZ_HIV_MTT , PYZ_pk, pkparam_MTT, FOCE())
#Inference
infer(pkfit_HIV_MTT )
#Calculate BIC
bic(pkfit_HIV_MTT)

#Testing HIV as a covariate on F
PYZ_HIV_F = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =1)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init = 45)
        """
        Mean transit time
        """
        tMTT ∈ RealDomain(; lower = 0, init = 2)
        """
        HIV F factor
        """
        HIV_F ∈ RealDomain(;  init = 0.5)
        """
          - ΩCL
          - ΩVc
          - ΩMTT
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [ 0.1, 0.1, 0.1, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.38)
    end
    @covariates begin
        """
        Weight (kg)
        """
        WT
        """
        HIV
        """
        HIV
    end
    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        CL = tvcl * (WT/median(PYZ_data.WT))^0.75 * exp(η[1]) 
        Vc = tvv * (WT/median(PYZ_data.WT)) *  exp(η[2])
        MTT = tMTT * exp(η[3])

        ktr = 3/MTT
        k20= CL/Vc
    end

    @dosecontrol begin
        bioav = ( ; Abs=1 * (1+ HIV_F*HIV) * exp(η[4])) 
    end

    @dynamics begin
        Abs' = -ktr*Abs
        B1' = ktr *Abs- ktr*B1
        B2' = ktr *B1- ktr*B2
        Central' = ktr*B2 - k20 * Central
    end
    @derived begin
        cp := @. log(0.00001+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p)
    end
end
#Initial parms
pkparam_F = init_params(PYZ_HIV_F )
#Fit using FOCE
pkfit_HIV_F = fit(PYZ_HIV_F , PYZ_pk, pkparam_F, FOCE())
#Inference
infer(pkfit_HIV_F)
#Calculate BIC
bic(pkfit_HIV_F)

############################################################################################
##########USING BEST FITTED 1 CMP AS THE BASE WITH WT #####################################
PYZ_cov = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =4)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init = 25)
        """
        Mean transit time
        """
        tMTT ∈ RealDomain(; lower = 0, init = 1)
        """
        Sex Factor CL, V, MTT and F (0 Male, 1 Female)
        """
        SexCL ∈ RealDomain(; init = 0.05)
        SexV ∈ RealDomain(; init = -0.2)
        SexMTT ∈ RealDomain(; init = 0.05)
        SexF ∈ RealDomain(;  init = 0.05)
        """
        AGE Factor CL, V, F and MTT 
        """
        AgeCL ∈ RealDomain(;  init = 0.05)
        AgeV ∈ RealDomain(;init = 0.05)
        AgeMTT ∈ RealDomain(;init = 0.05)
        AgeF ∈ RealDomain(;  init = 0.05)

        """
        Liver (ALTRES) CL
        """
        AltCL ∈ RealDomain(; init = 0.05)
        """
        KIDNEYS (CREATRES) CL
        """
        CreaCL ∈ RealDomain(; init = 0.05)
        """
        KIDNEYS (TOTBRES) CL
        """
        TotbCL ∈ RealDomain(; init = 0.05)
        """
          - ΩCL
          - ΩVc
          - ΩMTT
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [ 0.1, 0.1, 0.1, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.38)
    end
    @covariates begin
        """
        Weight (kg)
        """
        WT
        """
        Sex 
        """
        SEX
        """
        AGE
        """
        AGE
        """
        Liver function
        """
        ALTRES
        """
        Kidney function
        """
        CREATRES
        """
        Liver function
        """
        TOTBRES
    end
    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        CL = tvcl * (WT/median(PYZ_data.WT))^0.75 * (1+SexCL*SEX) * (1+AgeCL*(AGE-median(PYZ_data.AGE))) * (1+AltCL*(ALTRES-median(PYZ_data.ALTRES))) * (1+CreaCL*(CREATRES-median(PYZ_data.CREATRES))) * (1+TotbCL*(TOTBRES-median(PYZ_data.TOTBRES))) * exp(η[1]) 
        Vc = tvv * (WT/median(PYZ_data.WT)) *  (1+SexV*SEX) * (1+AgeV*(AGE-median(PYZ_data.AGE))) * exp(η[2])       
        MTT = tMTT * (1+SexMTT*SEX) * (1+AgeMTT*(AGE-median(PYZ_data.AGE))) * exp(η[3])

        ktr=3/MTT
        k20= CL/Vc
    end

    @dosecontrol begin
        bioav = ( ; Abs=1 * (1+SexF*SEX) * (1+AgeF*(AGE-median(PYZ_data.AGE))) * exp(η[4]))
    end

    @dynamics begin
        Abs' = -ktr*Abs
        B1' = ktr *Abs- ktr*B1
        B2' = ktr *B1- ktr*B2
        Central' = ktr*B2 - k20 * Central
    end
    @derived begin
        cp := @. log(0.00001+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p) 
    end
end

#Forward Covariate Selection
#######################################
#Initial parms
pkparam_for = init_params(PYZ_cov)
#Performing a forward covariate selection
covar_result = covariate_select(
    PYZ_cov,
    PYZ_pk,
    pkparam_for,
    FOCE();
    method = CovariateSelection.Forward,
    criterion = bic,
    control_param = (:SexCL, :SexV, :SexMTT, :SexF, :AgeCL, :AgeV, :AgeMTT, :AgeF, :AltCL, :CreaCL, :TotbCL) , # this is a Tuple
)
#Sorting results from best fit to worst
sort(DataFrame(covar_result.fits), :criterion)
#Printing the best model
covar_result.best_model
#Inference
infer(covar_result.best_model)
#Calculating BIC
bic(covar_result.best_model)
 
#BACKWARDS SELECTION, looking for 2 x log likelihood to increase by more than 6.635	(𝛼=0.01)
#Removing CREACL
PYZ_back1 =  @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =4)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init =25)
        """
        Mean transit time
        """
        tMTT ∈ RealDomain(; lower = 0, init = 1)
        """
          - ΩCL
          - ΩVc
          - ΩMTT
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [ 0.1, 0.1, 0.1, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.38)
    end
    @covariates begin
        """
        Weight (kg)
        """
        WT
    end
    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        CL = tvcl * (WT/median(PYZ_data.WT))^0.75 * exp(η[1]) 
        Vc = tvv * (WT/median(PYZ_data.WT)) * exp(η[2])
        MTT=tMTT * exp(η[3])

        ktr = 3/MTT 
        k20= CL/Vc
    end

    @dosecontrol begin
        bioav = ( ; Abs=1 * exp(η[4]))
    end

    @dynamics begin
        Abs' = -ktr*Abs
        B1' = ktr *Abs- ktr*B1
        B2' = ktr *B1- ktr*B2
        Central' = ktr*B2 - k20 * Central
    end
    @derived begin
        cp := @. log(0.0001+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p) 
    end
end
#Initial parms
pkparam_back1 = init_params(PYZ_back1 )
#Fit using FOCE
pkfit_back1= fit(PYZ_back1 , PYZ_pk, pkparam_back1, FOCE())
#Calculate BIC
bic(pkfit_back1)

#-2LOGLIK CALC
#Log likelihood reduced =-503.59135
#Log likelihood full= -489.18572
-2*(-503.59135-(-489.18572)) #28.811259999999947>6.635 , reduced model rejected

#Removing WT and adding back CreaCL
PYZ_back2 = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =4)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init =25)
        """
        Mean transit time
        """
        tMTT ∈ RealDomain(; lower = 0, init =1)
        """
        KIDNEYS (CREATRES) CL
        """
        CreaCL ∈ RealDomain(; init = 0.05)
        """
          - ΩCL
          - ΩVc
          - ΩMTT
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [ 0.1, 0.1, 0.1, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.38)
    end
    @covariates begin
        """
        Weight (kg)
        """
        WT
        """
        Kidney function
        """
        CREATRES
    end
    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        CL = tvcl * (1+CreaCL*(CREATRES-median(PYZ_data.CREATRES))) *  exp(η[1]) 
        Vc = tvv *  exp(η[2])
        MTT=tMTT * exp(η[3])

        ktr = 2/MTT 
        k20= CL/Vc
    end

    @dosecontrol begin
        bioav = ( ; Abs=1 * exp(η[4]))
    end

    @dynamics begin
        Abs' = -ktr*Abs
        B1' = ktr *Abs- ktr*B1
        Central' = ktr*B1 - k20 * Central
    end
    @derived begin
        cp := @. log(0.0001+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p)
    end
end
#Initial parms
pkparam_back2 = init_params(PYZ_back2 )
#Fit using FOCE
pkfit_back2= fit(PYZ_back2 , PYZ_pk, pkparam_back2, FOCE())
#Log likelihood reduced =-533.3601
#Log likelihood full model = -489.18572
-2*(-533.3601-(-489.18572)) #88.34875999999997>6.635 , reduced model rejected

#Both covariates included through back selection, CREATINE not clinically significant hence removed
#So through backwards selection best model only has WT as a covariate

#############################################################
#Now trying out TIME as a covariate (D0 compared to W6) on CL 
PYZ_TP = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =3)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init =35)
        """
        Mean transit time
        """
        tMTT ∈ RealDomain(; lower = 0, init = 1)
        """
        Time Clearance Factor
        """
        TPCL ∈ RealDomain(; init = 0.)
        """
          - ΩCL
          - ΩVc
          - ΩMTT
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [ 0.1, 0.1, 0.1, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.38)
    end
    @covariates begin
        """
        Weight (kg)
        """
        WT
        """
        Time (D0 =0, W6 =1)
        """
        TP
    end
    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        CL = tvcl * (WT/median(PYZ_data.WT))^0.75 * (1+TP*TPCL)* exp(η[1]) 
        Vc = tvv * (WT/median(PYZ_data.WT)) * exp(η[2])
        MTT=tMTT * exp(η[3])

        ktr = 3/MTT 
        k20= CL/Vc
    end

    @dosecontrol begin
        bioav = ( ; Abs=1 * exp(η[4]))
    end

    @dynamics begin
        Abs' = -ktr*Abs
        B1' = ktr *Abs- ktr*B1
        B2' = ktr *B1- ktr*B2
        Central' = ktr*B2 - k20 * Central
    end
    @derived begin
        cp := @. log(0.0001+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p) 
    end
end
#Initial parms
pkparam_TPCL = init_params(PYZ_TP )
#Fit using FOCE
pkfit_TP = fit(PYZ_TP , PYZ_pk, pkparam_TPCL, FOCE())
#Inference
infer(pkfit_TP)
#Calculating BIC
bic(pkfit_TP)
#Hence no improved fit

#############################################################
#Now trying out TIME as a covariate (D0 compared to W6) on F 
PYZ_TP_F = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =2)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init =50)
        """
        Mean transit time
        """
        tMTT ∈ RealDomain(; lower = 0, init = 1)
        """
        Time Bioavailability Factor
        """
        TPF ∈ RealDomain(; init = 0.)
        """
          - ΩCL
          - ΩVc
          - ΩMTT
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [ 0.1, 0.1, 0.1, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.38)
    end
    @covariates begin
        """
        Weight (kg)
        """
        WT
        """
        Time (D0 =0, W6 =1)
        """
        TP
    end
    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        CL = tvcl * (WT/median(PYZ_data.WT))^0.75 *  exp(η[1]) 
        Vc = tvv * (WT/median(PYZ_data.WT)) * exp(η[2])
        MTT=tMTT * exp(η[3])

        ktr = 3/MTT 
        k20= CL/Vc
    end

    @dosecontrol begin
        bioav = ( ; Abs=1 * (1+TP*TPF)* exp(η[4]))
    end

    @dynamics begin
        Abs' = -ktr*Abs
        B1' = ktr *Abs- ktr*B1
        B2' = ktr *B1- ktr*B2
        Central' = ktr*B2 - k20 * Central
    end
    @derived begin
        cp := @. log(0.0001+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p) #To get ng/mL
    end
end
#Initial parms
pkparam_TPF = init_params(PYZ_TP_F )
#Fit using FOCE
pkfit_TP_F = fit(PYZ_TP_F , PYZ_pk, pkparam_TPF, FOCE())
#Inference
infer(pkfit_TP_F)
#Calculate BIC
bic(pkfit_TP_F)
#Hence no improved fit

########################################################
#STEPWISE REMOVING OMEGAS - remove stepwise manually
PYZ_om1 = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =2)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init =50)
        """
        Mean transit time
        """
        tMTT ∈ RealDomain(; lower = 0, init = 1)
        """
          - ΩCL
          - ΩMTT
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.1, 0.1, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.18)
    end
    @covariates begin
        """
        Weight (kg)
        """
        WT
    end
    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        CL = tvcl * (WT/median(PYZ_data.WT))^0.75 * exp(η[1])
        Vc = tvv * (WT/median(PYZ_data.WT)) * exp(η[2])
        MTT=tMTT 

        ktr = 3/MTT * exp(η[3])
        k20= CL/Vc
    end

    @dosecontrol begin
        bioav = ( ; Abs=1)
    end

    @dynamics begin
        Abs' = -ktr*Abs
        B1' = ktr *Abs- ktr*B1
        B2' = ktr *B1- ktr*B2
        Central' = ktr*B2 - k20 * Central
    end
    @derived begin
        cp := @. log(0.0001+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p) #To get ng/mL
    end
end
#Initial Parms
pkparam_om1 = init_params(PYZ_om1 )
#Fit using FOCE
pkfit_PYZom1= fit(PYZ_om1 , PYZ_pk, pkparam_om1, FOCE())
#Inference
infer_PYZ=infer(pkfit_PYZom1)
#Calculate BIC
bic(pkfit_PYZom1)

#REMOVING OMV
#reduced model log-likelihood = -503.91634
#full model log-likelihood = -503.59135
-2*(-503.91634-(-503.59135)) #0.64998<6.635 , reduced model accepted

#REMOVING OMCL
#reduced model log-likelihood = -527.96522
#full model log-likelihood = -503.59135
-2*(-527.96522-(-503.59135)) #48.74774>6.635 , reduced model rejected

#REMOVING OMMTT
#reduced model log-likelihood = -650.66416
#full model log-likelihood = -503.59135
-2*(-650.66416-(-503.59135)) #294.14562>6.635 , reduced model rejected

#REMOVING OMF
#reduced model log-likelihood = -506.93106
#full model log-likelihood = -503.59135
-2*(-506.93106-(-503.59135)) #6.67942>6.635 , reduced model rejected

#Hence om on V removed, 
#Remove rest stepwise manually
PYZ_om2 = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =2)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init =48)
        """
        Mean transit time
        """
        tMTT ∈ RealDomain(; lower = 0, init = 1)
        """
          - ΩF
          - ΩMTT
        """
        Ω ∈ PDiagDomain(init = [0.1, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.18)
    end
    @covariates begin
        """
        Weight (kg)
        """
        WT
    end
    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        CL = tvcl * (WT/median(PYZ_data.WT))^0.75 * exp(η[1]) 
        Vc = tvv * (WT/median(PYZ_data.WT)) * exp(η[2])
        MTT=tMTT 

        ktr = 3/MTT 
        k20= CL/Vc
    end

    @dosecontrol begin
        bioav = ( ; Abs=1* exp(η[2]))
    end

    @dynamics begin
        Abs' = -ktr*Abs
        B1' = ktr *Abs- ktr*B1
        B2' = ktr *B1- ktr*B2
        Central' = ktr*B2 - k20 * Central
    end
    @derived begin
        cp := @. log(0.0001+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p) 
    end
end
#Initial parms
pkparam_om2 = init_params(PYZ_om2 )
#Fit using FOCE
pkfit_PYZom2= fit(PYZ_om2 , PYZ_pk, pkparam_om2, FOCE())
#Inference
infer_PYZ=infer(pkfit_PYZom2)
#Calculate BIC
bic(pkfit_PYZom2)

#REMOVING OMCL
#reduced model log-likelihood =  -600.13438
#full model log-likelihood (v removed) = -503.91634
-2*( -600.13438-(-503.91634)) #192.43608>6.635 , reduced model rejected

#REMOVING OMMTT
#reduced model log-likelihood = -654.64362
#full model log-likelihood = -503.91634
-2*(-654.64362-(-503.91634)) #301.45456>6.635 , reduced model rejected

#REMOVING OMF
#reduced model log-likelihood = -506.93106
#full model log-likelihood = -503.91634
-2*(-694.12091-(-503.91634)) #380.40914>6.635 , reduced model rejected

###############################################
#So the best model is:
PYZ = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =2)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init =50)
        """
        Mean transit time
        """
        tMTT ∈ RealDomain(; lower = 0, init = 1)
        """
          - ΩCL
          - ΩMTT
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.1, 0.1, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.18)
    end
    @covariates begin
        """
        Weight (kg)
        """
        WT
    end
    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        CL = tvcl * (WT/median(PYZ_data.WT))^0.75* exp(η[1])
        Vc = tvv * (WT/median(PYZ_data.WT)) 
        MTT=tMTT * exp(η[2])

        ktr = 3/MTT 
        k20= CL/Vc
    end

    @dosecontrol begin
        bioav = ( ; Abs=1 * exp(η[3]))
    end

    @dynamics begin
        Abs' = -ktr*Abs
        B1' = ktr *Abs- ktr*B1
        B2' = ktr *B1- ktr*B2
        Central' = ktr*B2 - k20 * Central
    end
    @derived begin
        cp := @. log(0.0001+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p) 
    end
end
#Initial parms
pkparam = init_params(PYZ )
#Fit using FOCE
pkfit_PYZ= fit(PYZ , PYZ_pk, pkparam, FOCE())
#Inference
infer_PYZ=infer(pkfit_PYZ)
#Calculate BIC
bic(pkfit_PYZ)
#Calculate shrinkage
ηshrinkage(pkfit_PYZ)
ϵshrinkage(pkfit_PYZ)
#Fitted parameter values
fittedparms =(tvcl=3.475,
              tvv = 34.623,
              tMTT = 0.56177, 
              Ω = [0.04 0.0 0.0; 0.0 0.26064 0.0; 0.0 0.0 0.015941],
              σ_p = 0.30388,)

#VPC for all data including months              
vpc_pk_final = vpc(
        pkfit_PYZ; 
        samples=200,
        covariates = [:tad],
        stratify_by = [:TP2]

        )

vpc_plot(
    PYZ,
    vpc_pk_final;
    figure = (; resolution = (1400, 1000), fontsize = 22),
    axis = (;
        xlabel = "Time after dose (hr)",
        ylabel = "Observed/Predicted\n Log Concentration (ng/mL)",
    ),
    facet = (; combinelabels = true),
)

#Removing month data for vpc of just D0 and W6
PYZ_D0W6 = filter(row -> row.TP2 in [0, 1], PYZ_data)
#New Pumas poulation for this data
PYZ_pk3 = read_pumas(
    PYZ_D0W6 ;
    id = :ID,
    time = :TIME,
    amt = :AMT,
    observations = [:LOG_DV],
    evid = :EVID,
    cmt = :CMT,
    covariates= [:WT, :HIV, :SEX, :ALTRES, :CREATRES, :AGE, :TOTBRES, :TP, :TP2]
)
#Refit to this data using FOCE
pkfit_PYZ2= fit(PYZ , PYZ_pk3, pkparam, FOCE())

#Produce vpc plots
vpc_pk_final = vpc(
        pkfit_PYZ2; 
        samples=2000,
        covariates = [:tad],
        stratify_by = [:TP2]

        )
vpc_plot(
    PYZ,
    vpc_pk_final;
    figure = (; resolution = (1400, 1000), fontsize = 22),
    axis = (;
        xlabel = "Time after dose (hr)",
        ylabel = "Observed/Predicted\n Log Concentration (ng/mL)",
    ),
    facet = (; combinelabels = true),
)

########## Create custom TAD GOF plot ##################################
#get inpsect element as dataframe
res_inspect_PYZ_df = DataFrame(res_inspect_PYZ)

#extract desired oclumns for custom plotting
gof_df = select(res_inspect_PYZ_df, [:LOG_DV_pred, :LOG_DV_ipred, :LOG_DV_wres, :LOG_DV_iwres, :tad])

# Remove missing values
gof_df = dropmissing(gof_df)

# Map to Float32 for scatter plot
gof_df = mapcols(x -> Float32.(x), gof_df)

# OLS requires Float64 type
df_ols = mapcols(x -> Float64.(x), gof_df)

# Create OLS model
tad_wres_ols_model = lm(@formula(LOG_DV_wres ~ tad), df_ols)
tad_wres_ols_fit = predict(tad_wres_ols_model)

# Create a new df for Loess to ensure data is ordered by :tad
df_loess_sorted = sort(gof_df, :tad)

# Create Loess model
tad_wres_loess_model = loess(df_loess_sorted.tad, df_loess_sorted.LOG_DV_wres)

# Generate predictions using the Loess model
x_vals = LinRange(minimum(df_loess_sorted.tad), maximum(df_loess_sorted.tad), 100)
tad_wres_loess_fit = Loess.predict(tad_wres_loess_model, x_vals)

# Create figure for plot
tad_fig = Figure(resolution = (800, 600), title = "LNDV")

# Create an axis
ax = Axis(tad_fig[1, 1], 
    title = "LOG_DV\nConcentration (ng/ml)", 
    titlesize = 20, 
    titlegap = 10, 
    xlabel = "Time after dose", 
    ylabel = "Weighted population residuals"
)

# Add scatter points of WRES and TAD
CairoMakie.scatter!(ax, gof_df.tad, gof_df.LOG_DV_wres, color = :black)

# Add OLS line
lines!(ax, df_ols.tad, tad_wres_ols_fit, color = :green, linewidth = 2)

# Add Loess line
lines!(ax, x_vals, tad_wres_loess_fit, color = :red, linewidth = 2)

# Display the figure
display(tad_fig)

#Goodness-of-fit Plot
res_inspect_PYZ= inspect(pkfit_PYZ)
gof_PYZ = goodness_of_fit(res_inspect_PYZ; figure = (; fontsize = 12))
"julia.execution.resultType": "inline",
"julia.plotpane.show": true
#Checking for correlation
correlation_diagnostic(
    infer_PYZ)

#500 sample non-parametric bootstrap
boot= infer(
    pkfit_PYZ,
    Bootstrap(
        samples = 500,
        ensemblealg = EnsembleThreads(),
    ),
)
#Load results into a dataframe
boot_PYZ=DataFrame(boot.vcov)

#Produce density plots of parameter estimates
StatsPlots.density(boot_PYZ.tvcl, label=false, xlabel="tvcl", ylabel="Density", title="Density Plot of tvcl", Legend=false)
StatsPlots.density(boot_PYZ.tvv, label=false, xlabel="tvv", ylabel="Density", title="Density Plot of tvv", Legend=false)
StatsPlots.density(boot_PYZ.tMTT, label=false, xlabel="tMTT", ylabel="Density", title="Density Plot of tMTT", Legend=false)
StatsPlots.density(boot_PYZ."Ω₃,₃", label=false, xlabel="Ω₃,₃", ylabel="Density", title="Density Plot of Ω₃,₃", Legend=false)

#Summary stats of bootstrap
summary_stats = describe(boot_PYZ)

####################
#Calculating skewness and kurtosis
#Initialise vectors to calculate skewness and kurtosis
parameters = String[]
skewness_vals = Float64[]
kurtosis_vals = Float64[]
# Calculate skewness and kurtosis for each parameter
for col in names(boot_PYZ)
    push!(parameters, col)
    push!(skewness_vals, skewness(boot_PYZ[!, col]))
    push!(kurtosis_vals, kurtosis(boot_PYZ[!, col]))
end
#Produce dataframe of results
skewness_kurtosis = DataFrame(Parameter = parameters, Skewness = skewness_vals, Kurtosis = kurtosis_vals)

######################
#Calculating Bootstrap Bias and Variance
original_estimates = Dict("tvcl" => 2, "tvv" => 50, "tMTT" => 1, "Ω₁,₁" => 0.1, "Ω₂,₂" => 0.1, "Ω₃,₃" => 0.1, "σ_p" => 0.18)
# Initialize vectors to store bias and variance
parameters = String[]
bias_vals = Float64[]
variance_vals = Float64[]
# Calculate bias and variance for each parameter
for col in names(boot_PYZ)
    original_value = original_estimates[col]
    bootstrap_mean = mean(boot_PYZ[!, col])
    bias = bootstrap_mean - original_value
    variance = var(boot_PYZ[!, col])
    push!(parameters, col)
    push!(bias_vals, bias)
    push!(variance_vals, variance)
end
# Create DataFrame
bootstrap_bias_variance = DataFrame(Parameter = parameters, Bias = bias_vals, Variance = variance_vals)

################################
#Conduction Kolmogorov-Smirnov test, testing for normality
# Initialize vectors to store Kolmogorov-Smirnov test results
for col in names(boot_PYZ)
    boot_PYZ[!, col] = convert(Vector{Float64}, boot_PYZ[!, col])
end
parameters = String[]
ks_p_vals = Float64[]

# Perform Kolmogorov-Smirnov test
for col in names(boot_PYZ)
    test = HypothesisTests.ApproximateOneSampleKSTest(boot_PYZ[!, col], Normal(mean(boot_PYZ[!, col]), std(boot_PYZ[!, col])))
    push!(parameters, col)
    push!(ks_p_vals, HypothesisTests.pvalue(test))
end

# Create DataFrame for Kolmogorov-Smirnov test results
ks_results = DataFrame(Parameter = parameters, p_value = ks_p_vals)

# Generate Q-Q plots
for col in names(boot_PYZ)
    plot = @df boot_PYZ StatsPlots.qqplot(Normal(mean(boot_PYZ[!, col]), std(boot_PYZ[!, col])), boot_PYZ[!, col], title = "$col Q-Q Plot", legend = false)
    display(plot)
end

################################################
#DIAGNOSTIC PLOTS FOR D0
#D0 fit using FOCE
pkfit_PYZD0 = fit(PYZ, PYZ_pkD0, pkparam, FOCE())
#Inference
infer_PYZD0=infer(pkfit_PYZD0)
#Calculate NIC
bic(pkfit_PYZD0)
#Sims, VPC and Gof Plots
simpk_iparams_PYZDO = simobs(PYZ, PYZ_pkD0, fittedparms)
sim_plot(
    PYZ,
    simpk_iparams_PYZDO;
    observations = [:LOG_DV],
    figure = (; fontsize = 18),
    axis = (;
        xlabel = "Time (hr)",
        ylabel = "Observed/Predicted \n Log Concentration (ng/mL)",
    ),
)
pk_vpc_PYZD0 = vpc(
    pkfit_PYZD0,
    200;
    observations = [:LOG_DV],
    ensemblealg = EnsembleThreads(), # multi-threading
)
vpc_plot(
    PYZ,
    pk_vpc_PYZD0;
    figure = (; resolution = (1400, 1000), fontsize = 22),
    axis = (;
        xlabel = "Time (hr)",
        ylabel = "Observed/Predicted\n Log Concentration (ng/mL)",
    ),
    facet = (; combinelabels = true),
)

res_inspect_PYZD0 = inspect(pkfit_PYZD0)
gof_PYZD0 = goodness_of_fit(res_inspect_PYZD0; figure = (; fontsize = 12))

##############################################
################################################
#DIAGNOSTIC PLOTS FOR W6
#################################################
#Fit model using FOCE
pkfit_PYZW6 = fit(PYZ, PYZ_pkW6, pkparam, FOCE())
#Inference
infer_PYZW6 = infer(pkfit_PYZW6)
#Calculate BIC
bic(pkfit_PYZW6)
#VPC, Sims and GOF plots
simpk_iparams_PYZW6 = simobs(PYZ, PYZ_pkW6, fittedparms)
sim_plot(
    PYZ,
    simpk_iparams_PYZW6;
    observations = [:LOG_DV],
    figure = (; fontsize = 18),
    axis = (;
        xlabel = "Time (hr)",
        ylabel = "Observed/Predicted \n Log Concentration (ng/mL)",
    ),
)
pk_vpc_PYZW6 = vpc(
    pkfit_PYZW6,
    200;
    observations = [:LOG_DV],
    ensemblealg = EnsembleThreads(), # multi-threading
)
vpc_plot(
    PYZ ,
    pk_vpc_PYZW6;
    figure = (; resolution = (1400, 1000), fontsize = 22),
    axis = (;
        xlabel = "Time (hr)",
        ylabel = "Observed/Predicted\n Log Concentration (ng/mL)",
    ),
    facet = (; combinelabels = true),
)
res_inspect_PYZW6 = inspect(pkfit_PYZW6)
gof_PYZW6 = goodness_of_fit(res_inspect_PYZW6; figure = (; fontsize = 12))

####################################################
###############RUNNING INDIVIDUAL PLOTS#############
####################################################
pkparam = (; #initial parameters as model parameter estimates
tvcl = 3.475, tvv = 34.623, tMTT = 0.56177, 
Ω = Diagonal([0.000000001, 0.000000001, 0.000000001, 0.000000001]), #using zeroed etas
σ_p = 0.0000001)
using Makie
#Weight 21kg to 34kg
sim1 = DosageRegimen(800000000, cmt = 1, time = 0, ii = 24, addl = 100)
#Dosage regime of 100 800mg doses of pyrazinamide given every 24 hours
#adjsuting proportional and random error to 0 for visualisation purposes
patient21 = Subject(id = 60, 
                        events = sim1,
                        covariates = (WT = 21,),  
                        observations = (LNDV = nothing,)
                        )
patient34 = Subject(id = 61,
                        events = sim1,
                        covariates = (WT = 34,),  
                        observations = (LNDV = nothing,)
                        )
Random.seed!(123)
sim21 = simobs(PYZ, patient21, pkparam, obstimes = 0:1:1008)
sim34 = simobs(PYZ, patient34, pkparam, obstimes = 0:1:1008)

# Create a new figure
fig = Figure()
# Add an axis to the figure
ax = Axis(fig[1, 1], title = "Simulation Plot", xlabel = "Time", ylabel = "LOG_DV")
# Plot the first simulation data
lines!(ax, sim21.time, sim21.observations[:LOG_DV]; color = :mediumpurple, label = "Simulation of patient weight 21kg")
# Plot the second simulation data
lines!(ax, sim34.time, sim34.observations[:LOG_DV]; color = RGB(1.0, 0.6, 0.4), label = "Simulation of patient weight 34kg")
# Add legend
axislegend(ax, position=(1, 0))  # (1, 0)
# Display the figure
fig

#Weight 35kg to 39kg
sim1 = DosageRegimen(1000000000, cmt = 1, time = 0, ii = 24, addl = 100)
#Dosage regime of 100 1000mg doses of pyrazinamide given every 24 hours
#adjsuting proportional and random error to 0 for visualisation purposes
patient35 = Subject(id = 60,
                        events = sim1,
                        covariates = (WT = 35,),  
                        observations = (LNDV = nothing,)
                        )
patient39 = Subject(id = 61,
                        events = sim1,
                        covariates = (WT = 39,),  
                        observations = (LNDV = nothing,)
                        )
Random.seed!(123)
sim35 = simobs(PYZ, patient35, pkparam, obstimes = 0:1:1008)
sim39 = simobs(PYZ, patient39, pkparam, obstimes = 0:1:1008)

# Create a new figure
fig = Figure()
# Add an axis to the figure
ax = Axis(fig[1, 1], title = "Simulation Plot", xlabel = "Time", ylabel = "LOG_DV")
# Plot the first simulation data
lines!(ax, sim35.time, sim35.observations[:LOG_DV]; color = :mediumpurple, label = "Simulation of patient weight 35kg")
# Plot the second simulation data
lines!(ax, sim39.time, sim39.observations[:LOG_DV]; color = RGB(1.0, 0.6, 0.4), label = "Simulation of patient weight 39kg")
# Add legend
axislegend(ax, position=(1, 0))  # (1, 0)
# Display the figure
fig

#Weight 40kg to 54kg
sim1 = DosageRegimen(1200000000, cmt = 1, time = 0, ii = 24, addl = 100)
#Dosage regime of 100 1200mg doses of pyrazinamide given every 24 hours
#adjsuting proportional and random error to 0 for visualisation purposes
patient40 = Subject(id = 60,
                        events = sim1,
                        covariates = (WT = 35,),  
                        observations = (LNDV = nothing,)
                        )
patient54 = Subject(id = 61,
                        events = sim1,
                        covariates = (WT = 39,),  
                        observations = (LNDV = nothing,)
                        )
Random.seed!(123)
sim40 = simobs(PYZ, patient40, pkparam, obstimes = 0:1:1008)
sim54 = simobs(PYZ, patient54, pkparam, obstimes = 0:1:1008)

# Create a new figure
fig = Figure()
# Add an axis to the figure
ax = Axis(fig[1, 1], title = "Simulation Plot", xlabel = "Time", ylabel = "LOG_DV")
# Plot the first simulation data
lines!(ax, sim40.time, sim40.observations[:LOG_DV]; color = :mediumpurple, label = "Simulation of patient weight 40kg")
# Plot the second simulation data
lines!(ax, sim54.time, sim54.observations[:LOG_DV]; color = RGB(1.0, 0.6, 0.4), label = "Simulation of patient weight 54kg")
# Add legend
axislegend(ax, position=(1, 0))  # (1, 0)
# Display the figure
fig

#Weight 55kg to 70kg
sim1 = DosageRegimen(1600000000, cmt = 1, time = 0, ii = 24, addl = 1)
#Dosage regime of 100 1600mg doses of pyrazinamide given every 24 hours
#adjsuting proportional and random error to 0 for visualisation purposes
patient55 = Subject(id = 60,
                        events = sim1,
                        covariates = (WT = 35,),  
                        observations = (LNDV = nothing,)
                        )
patient70 = Subject(id = 61,
                        events = sim1,
                        covariates = (WT = 39,),  
                        observations = (LNDV = nothing,)
                        )
Random.seed!(123)
sim55 = simobs(PYZ, patient55, pkparam, obstimes = 0:1:24)
sim70 = simobs(PYZ, patient70, pkparam, obstimes = 0:1:24)
maximum(PYZ_data.WT)
PYZ_data.WT
# Create a new figure
fig = Figure()
# Add an axis to the figure
ax = Axis(fig[1, 1], title = "Simulation Plot", xlabel = "Time", ylabel = "LOG_DV")
# Plot the first simulation data
lines!(ax, sim55.time, sim55.observations[:LOG_DV]; color = :mediumpurple, label = "Simulation of patient weight 55kg")
# Plot the second simulation data
lines!(ax, sim70.time, sim70.observations[:LOG_DV]; color = RGB(1.0, 0.6, 0.4), label = "Simulation of patient weight 70kg")
# Add legend
axislegend(ax, position=(1, 0))  # (1, 0)
# Display the figure
fig

#Calculate predicted CMAX from each simulated weight group
cmax_21 = maximum(sim21.observations[:LOG_DV])
cmax_34 = maximum(sim34.observations[:LOG_DV])
cmax_35 = maximum(sim35.observations[:LOG_DV])
cmax_39 = maximum(sim39.observations[:LOG_DV])
cmax_40 = maximum(sim40.observations[:LOG_DV])
cmax_54 = maximum(sim54.observations[:LOG_DV])
cmax_55 = maximum(sim55.observations[:LOG_DV])
cmax_70 = maximum(sim70.observations[:LOG_DV])


#CALCULATING AUC PREDICTIONS
weight_bands = [
    (21, 34, 800),
    (35, 39, 1000),
    (40, 54, 1200),
    (55, 70, 1600),
    (71, 85, 2000) # Assuming an upper limit of 100 kg for simplicity
]

# Initialize a structure to hold AUC results
auc_results_WT = OrderedDict{String, Vector{Float64}}()

# Function to calculate AUC
function calculate_auc(sim_data)
    auc = 0.0
    time_points = sim_data.time
    concentrations = sim_data.observations.LOG_DV
    for i in 1:(length(time_points) - 1)
        dt = time_points[i+1] - time_points[i]
        auc += (exp(concentrations[i]) + exp(concentrations[i+1])) / (2*1000) * dt #mg/L
    end
    return auc
end

# Run simulations for each weight band
for (lower, upper, dose) in weight_bands
    weight_band_label = "$lower-$upper kg"
    auc_results_WT[weight_band_label] = Float64[]
    
    for _ in 1:100
        # Generate a random weight within the band
        weight = rand(lower:upper)
        
        # Define the dosage regimen
        regimen = DosageRegimen(dose * 1e6, cmt = 1, time = 0, ii = 24, addl = 1)
        
        # Create a subject with the given weight
        subject = Subject(id=1, events=regimen, covariates=(WT=weight,))
        
        # Simulate the model
        sim = simobs(PYZ, subject, init_params(PYZ), obstimes=0:1:24)
        
        # Calculate AUC for the simulated data
        auc = calculate_auc(sim)

        push!(auc_results_WT[weight_band_label], auc)
    end
end

#Get keys in the weight band
weight_bands_labels = collect(keys(auc_results_WT))
#Assign to data for plotting
auc_data_WT = [auc_results_WT[band] for band in weight_bands_labels]

# Create the grouped box plot with weight bands on the x-axis
Plots.boxplot(
    auc_data_WT,             # Data for boxplots (AUC values)
    xlabel = "Weight Band (kg)",
    ylabel = "AUC (mg/L)",
    title = "AUC for Different Weight Bands",
    box=:outliers,        # Show outliers
    orientation = :vertical,  # Weight bands on x-axis
    xticks = (1:length(weight_bands_labels), weight_bands_labels),  # Set x-axis ticks and labels
    legend = false        # No legend
)

#Print summary stats
describe.(auc_data_WT)

#Produce density plots
density_plots = Plots.plot(title="Density Plot of AUC for Different Weight Bands")
for (label, auc_values) in auc_results_WT
    Plots.density!(auc_values, label=label)
end
density_plots

#CDF Plot
function plot_cdf(data, labels, title)
    p = Plots.plot(title=title, xlabel="AUC (mg/L)", ylabel="CDF", legend=:bottomright)
    for (label, auc_values) in zip(labels, data)
        ecdf_func = ecdf(auc_values)
        x_values = sort(auc_values)
        y_values = ecdf_func.(x_values)
        Plots.plot!(p, x_values, y_values, label=label)
    end
    return p
end
# Plot CDF
cdf_plot = plot_cdf(auc_data_WT, weight_bands_labels, "CDF of AUC for Different Weight Bands")


#Violin Plots
# Flatten the data and create a corresponding label array
flattened_data = vcat(auc_data_WT...)
labels = repeat(weight_bands_labels, inner=100)
# Generate the violin plot
StatsPlots.violin(labels, flattened_data, xlabel="Weight Band (kg)", ylabel="AUC (mg/L)", title="Violin Plot of AUC for Different Weight Bands", legend=false)

# Prepare data for Kruskal-Wallis test
data = vcat([auc_results_WT[band] for band in weight_bands_labels]...)
labels = repeat(weight_bands_labels, inner=100)

# Perform the Kruskal-Wallis test
groups = unique(labels)
grouped_data = [flattened_data[labels .== group] for group in groups]
kw_test = HypothesisTests.KruskalWallisTest(grouped_data...)

################################
#AUC for observed data by weight
weight_bands = [
    (21, 34, "21-34 kg"),
    (35, 39, "35-39 kg"),
    (40, 54, "40-54 kg"),
    (55, 70, "55-70 kg"),
]

#Function to calculate AUC for observed data on D0 and W6 then find the average
function calculate_auc_obs(time_points, dv)
    auc = 0.0
    for i in 1:(length(time_points) - 1)
        dt = time_points[i+1] - time_points[i]
        auc += ((dv[i]) + (dv[i+1])) / (2*1000) * dt #mg/L
    end
    return auc
end

# Initialize a dictionary to store AUC results for each weight band
auc_results_WT_obs = OrderedDict{String, Vector{Float64}}()

# Initialize AUC calculation for each weight band
for (lower, upper, label) in weight_bands
    auc_results_WT_obs[label] = Float64[]
end

# Remove rows with missing values
cleaned_data = dropmissing(PYZ_data)
# Group data by subjects and calculate AUC for each subject
grouped_data = groupby(cleaned_data, [:ID])

#Apply AUC function to each ID
for (id, group) in pairs(grouped_data)
    weight = group.WT[1]
    auc_tp2_0 = nothing
    auc_tp2_1 = nothing

    # Check for AUC with TP2=0
    sub_group_0 = filter(row -> row[:TP2] == 0, group)
    if !isempty(sub_group_0)
        auc_tp2_0 = calculate_auc_obs(sub_group_0.TAD, sub_group_0.DV)
    end

    # Check for AUC with TP2=1
    sub_group_1 = filter(row -> row[:TP2] == 1, group)
    if !isempty(sub_group_1)
        auc_tp2_1 = calculate_auc_obs(sub_group_1.TAD, sub_group_1.DV)
    end

    # Calculate the average AUC if any AUC values are present
    auc_values = filter(!isnothing, [auc_tp2_0, auc_tp2_1])
    if !isempty(auc_values)
        average_auc = mean(auc_values)
        for (lower, upper, label) in weight_bands
            if lower <= weight <= upper
                push!(auc_results_WT_obs[label], average_auc)
            end
        end
    end
end
#Collect weight bands label
weight_bands_labels = collect(keys(auc_results_WT_obs))
#Push to a dataframe for analysis
auc_data_WT_obs = [auc_results_WT_obs[band] for band in weight_bands_labels]

#Print summary stats
describe.(auc_data_WT_obs)

#HIV AUC##############################
hiv_bands = [
    (0, "HIV Negative"),
    (1, "HIV Positive")
]

# Initialize a dictionary to store AUC results for each HIV status band
auc_results_HIV_obs = OrderedDict{String, Vector{Float64}}()

# Initialize AUC calculation for each HIV status band
for (status, label) in hiv_bands
    auc_results_HIV_obs[label] = Float64[]
end

# Group data by subjects and calculate AUC for each subject
for (id, group) in pairs(grouped_data)
    hiv_status = group.HIV[1]
    auc_tp2_0 = nothing
    auc_tp2_1 = nothing

    # Check for AUC with TP2=0
    sub_group_0 = filter(row -> row[:TP2] == 0, group)
    if !isempty(sub_group_0)
        auc_tp2_0 = calculate_auc_obs(sub_group_0.TAD, sub_group_0.DV)
    end

    # Check for AUC with TP2=1
    sub_group_1 = filter(row -> row[:TP2] == 1, group)
    if !isempty(sub_group_1)
        auc_tp2_1 = calculate_auc_obs(sub_group_1.TAD, sub_group_1.DV)
    end

    # Calculate the average AUC if any AUC values are present
    auc_values = filter(!isnothing, [auc_tp2_0, auc_tp2_1])
    if !isempty(auc_values)
        average_auc = mean(auc_values)
        for (status, label) in hiv_bands
            if hiv_status == status
                push!(auc_results_HIV_obs[label], average_auc)
            end
        end
    end
end

#Collect HIV labels
hiv_bands_labels = collect(keys(auc_results_HIV_obs))
#Push results to a dataframe for analysis
auc_data_HIV_obs = [auc_results_HIV_obs[band] for band in hiv_bands_labels]

#Print summary stats
describe.(auc_data_HIV_obs)

#Produce a density plot stratified by HIV status
density_plots = Plots.plot(title="Density Plot of AUC for HIV Status")
for (label, auc_values) in auc_results_HIV_obs
    Plots.density!(auc_values, label=label)
end
density_plots

# Create the grouped box plot wstratified by HIV status
Plots.boxplot(
    auc_data_HIV_obs,             # Data for boxplots (AUC values)
    xlabel = "HIV Status",
    ylabel = "AUC (mg/L)",
    title = "AUC for HIV",
    box=:outliers,        # Show outliers
    orientation = :vertical,  # HIV status on x-axis
    xticks = (1:length(hiv_bands_labels), hiv_bands_labels),  # Set x-axis ticks and labels
    legend = false        # No legend
)

#CDF Plot
# Plot CDF
cdf_plot = plot_cdf(auc_data_HIV_obs, hiv_bands_labels, "CDF of AUC for HIV Status")

# Prepare data for Kruskal-Wallis test
flattened_data = vcat(auc_data_HIV_obs...)
labels = []

for (label, aucs) in zip(hiv_bands_labels, auc_data_HIV_obs)
    append!(labels, repeat([label], inner=length(aucs)))
end

# Ensure the lengths match
@assert length(flattened_data) == length(labels)

# Group data by labels
grouped_data = Dict{String, Vector{Float64}}()
for label in hiv_bands_labels
    grouped_data[label] = []
end
for (data, label) in zip(flattened_data, labels)
    push!(grouped_data[label], data)
end

# Prepare data for the Kruskal-Wallis test
kw_data = [grouped_data[label] for label in hiv_bands_labels]

# Perform the Kruskal-Wallis test
kw_test = KruskalWallisTest(kw_data...)

# Generate the violin plot stratified by HIV status
StatsPlots.violin(labels, flattened_data, xlabel="HIV Status(kg)", ylabel="AUC (mg/L)", title="Violin Plot of AUC for HIV Status", legend=false)

#Comparing D0 and W6 AUC##################
##########################################
tp2_bands = [(0, "D0"), (1, "W6")] 

# Initialize a dictionary to store AUC results for each TP2 status band
auc_results_TP2 = OrderedDict{String, Vector{Float64}}()

# Initialize AUC calculation for each TP2 status band
for (status, label) in tp2_bands
    auc_results_TP2[label] = Float64[]
end

# Group data by subjects and calculate AUC for each subject on D0 then again on W6
grouped_data = groupby(cleaned_data, [:ID, :TP2])

#Apply the function to each ID
for (key, group) in pairs(grouped_data)
    id, tp2_status = key
    if tp2_status in [0, 1]
        for (status, label) in tp2_bands
            if tp2_status == status
                # Calculate AUC for the subject
                auc = calculate_auc(group.TAD, group.DV)
                
                # Store the AUC in the corresponding TP2 status band
                push!(auc_results_TP2[label], auc)
            end
        end
    end
end
# Prepare data for analysis
tp2_bands_labels = collect(keys(auc_results_TP2))
auc_data_TP2 = [auc_results_TP2[band] for band in tp2_bands_labels]

#Print summary stats
describe.(auc_data_TP2)

# Density plots for AUC for TP2 status
density_plots = Plots.plot(title="Density Plot of AUC for Time Point")
for (label, auc_values) in auc_results_TP2
    Plots.density!(auc_values, label=label)
end
density_plots

# Create the grouped box plot with TP2 status on the x-axis
Plots.boxplot(
    auc_data_TP2,             # Data for boxplots (AUC values)
    xlabel = "Time Point Status",
    ylabel = "AUC (mg/L)",
    title = "AUC for Different Time Point Status",
    box=:outliers,        # Show outliers
    orientation = :vertical,  # TP2 status on x-axis
    xticks = (1:length(tp2_bands_labels), tp2_bands_labels),  # Set x-axis ticks and labels
    legend = false        # No legend
)

# CDF Plot
cdf_plot = plot_cdf(auc_data_TP2, tp2_bands_labels, "CDF of AUC for Different Time Point Status")

# Prepare data for Kruskal-Wallis test
flattened_data = vcat(auc_data_TP2...)
labels = []
for (label, aucs) in zip(tp2_bands_labels, auc_data_TP2)
    append!(labels, repeat([label], inner=length(aucs)))
end

# Ensure the lengths match
@assert length(flattened_data) == length(labels)

# Group data by labels
grouped_data = Dict{String, Vector{Float64}}()
for label in tp2_bands_labels
    grouped_data[label] = []
end
for (data, label) in zip(flattened_data, labels)
    push!(grouped_data[label], data)
end

# Prepare data for the Kruskal-Wallis test
kw_data = [grouped_data[label] for label in tp2_bands_labels]

# Perform the Kruskal-Wallis test
kw_test = KruskalWallisTest(kw_data...)

# Generate the violin plot
StatsPlots.violin(labels, flattened_data, xlabel="Time Point Status", ylabel="AUC (mg/L)", title="Violin Plot of AUC for Different Time Point Status", legend=false)

#AUC results for whole data, no stratification#
# Initialize dictionary to store AUC results
auc_results_obs = Float64[]
# Calculate AUC for each ID
grouped_data = groupby(cleaned_data, [:ID, :TP2])
for (key, group) in pairs(grouped_data)
    id, tp2_status = key
    if tp2_status in [0, 1]
        # Calculate AUC for the subject
        auc = calculate_auc(group.TAD, group.DV)
        
        # Store the AUC in the results vector
        push!(auc_results_obs, auc)
    end
end
#Print summary stats
describe(auc_results_obs)
