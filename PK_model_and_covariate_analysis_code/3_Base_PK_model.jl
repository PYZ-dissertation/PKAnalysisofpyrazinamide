############
#Julia code to test for the best fitted PK model to the pyrazinamide concentration data
#using bayesian inference criterion 

#Models tested:
# ~ One-compartment first order absorption
# ~ One-compartment with 1-4 transit compartments (4th for completeness)
# ~ Two-compartment with 1-2 transit compartments

#Definitions:
# ~ CL - central compartment clearance (L/h)
# ~ CLp - peripheral compartment clearnace (L/h)
# ~ Vc - central compartment apparant volume of distribution (L)
# ~ Vp - peripheral compartment apparant volume of distribution (L)
# ~ MTT - mean transit time (h^(-1))
# ~ ka - first-order absorption rate (h^(-1))
# ~ ktr - transit compartment rate ((n+1)/MTT where n is the number of compartments, not including absorption compartments)
# ~ η - IIV of each parameter
############

#Load packages
using XLSX, Pumas, Dates, DataFrames, DataFramesMeta, DataStructures

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
    if col in [ "TOTBRES", "CREATRES", "WT"]
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

# Adding :DATE and :TIMEC to create :CTIM
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

filter!((row) -> row.TAD ≤ 35, PYZ_data) #ID 27 has some values misspecified, removing those

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

# Apply the function to each group and update the TIME column
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
# Add TP column to PYZ_data
PYZ_data = add_tp_column!(PYZ_data)

# Set ROUTE column to "ev" for all rows
PYZ_data.ROUTE .= "ev"

#Data now in correct format, can begin PK analysis
#####################

#Read in a Pumas model for the data
PYZ_pk = read_pumas(
    PYZ_data;
    id = :ID,
    time = :TIME,
    amt = :AMT,
    observations = [:LOG_DV],
    evid = :EVID,
    cmt = :CMT,
    covariates= [:WT, :HIV, :SEX, :ALTRES, :CREATRES, :AGE, :TOTBRES, :TP],
    route =:ROUTE,
)

###########PK MODELLING#########################
#ONE-COMPARTMENT MODEL, 1ST ORDER ABSORPTION####
pk_1cmp_abs = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =4)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init = 35)
        """
        Absorption rate
        """
        tka ∈ RealDomain(; lower = 0, init = 1)
        """
          - ΩCL
          - ΩVc
          - Ωka
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [ 0.1, 0.1, 0.1, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.38)
    end

    @random begin
        η ~ MvNormal(Ω)
    end
    @pre begin
        CL = tvcl * exp(η[1])
        Vc = tvv * exp(η[2])

        ka = (tka) * exp(η[3])
        k20= CL/Vc
    end
    @dosecontrol begin
        bioav = ( ; Abs=1 * exp(η[4]))
    end

    @dynamics begin
        Abs' = -ka*Abs
        Central' = ka*Abs - k20 * Central
    end
    @derived begin
        cp := @. log(0.01+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p) #To get ng/mL
    end

end
#Initial parms
pkparam_1abs = init_params(pk_1cmp_abs)
#Fitting using FOCE
pkfit_1cmp_abs = fit(pk_1cmp_abs, PYZ_pk, pkparam_1abs, FOCE())
#Inference
infer(pkfit_1cmp_abs)
#Calc BIC
bic(pkfit_1cmp_abs)

###ONE-COMPARTMENT 1 TRANS####
pk_1cmp_tran1 = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =3)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init = 30)
        """
        Mean transit time
        """
        tMTT ∈ RealDomain(; lower = 0, init = 2)
        """
          - ΩCL
          - ΩVc
          - ΩMTT
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [ 0.1, 0.1, 0.1,0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.38)
    end

    @random begin
        η ~ MvNormal(Ω)
    end
    @pre begin
        CL = tvcl * exp(η[1])
        Vc = tvv * exp(η[2])
        MTT = tMTT * exp(η[3])

        ktr = 2/ MTT
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
        cp := @. log(0.1+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p) 
    end

end
#Initial parms
pkparam_1tran = init_params(pk_1cmp_tran1)
#Fit using FOCE
pkfit_1cmp_tran1 = fit(pk_1cmp_tran1, PYZ_pk, pkparam_1tran, FOCE())
#Inference
infer(pkfit_1cmp_tran1)
#Calculate BIC
bic(pkfit_1cmp_tran1)

###ONE-COMPARTMENT 2 TRANS####
pk_1cmp_tran2 = @model begin
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

    @random begin
        η ~ MvNormal(Ω)
    end
    @pre begin
        CL = tvcl * exp(η[1])
        Vc = tvv * exp(η[2])
        MTT = tMTT * exp(η[3])

        ktr = 3/ MTT
        k20= CL/Vc
    end
    @dosecontrol begin
        bioav = ( ; Abs=1 * exp(η[4]))
    end
    @dynamics begin
        Abs' = -ktr*Abs
        B1' = ktr *Abs- ktr*B1
        B2' = ktr* B1 - ktr*B2
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
pkparam_1tran2 = init_params(pk_1cmp_tran2)
#Fit using FOCE
pkfit_1cmp_tran2 = fit(pk_1cmp_tran2, PYZ_pk, pkparam_1tran2, FOCE())
#Inference
infer(pkfit_1cmp_tran2)
#Calculate BIC
bic(pkfit_1cmp_tran2)

###ONE-COMPARTMENT 3 TRANS####
pk_1cmp_tran3 = @model begin
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

    @random begin
        η ~ MvNormal(Ω)
    end
    @pre begin
        CL = tvcl * exp(η[1])
        Vc = tvv * exp(η[2])
        MTT = tMTT * exp(η[3])

        ktr = 4/ MTT
        k20= CL/Vc
    end
    @dosecontrol begin
        bioav = ( ; Abs=1 * exp(η[4]))
    end
    @dynamics begin
        Abs' = -ktr*Abs
        B1' = ktr *Abs- ktr*B1
        B2' = ktr* B1 - ktr*B2
        B3' = ktr* B2 -ktr*B3
        Central' = ktr*B3 - k20 * Central
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
pkparam_1tran3 = init_params(pk_1cmp_tran3)
#Fit using FOCE
pkfit_1cmp_tran3 = fit(pk_1cmp_tran3, PYZ_pk, pkparam_1tran3, FOCE())
#Inference
infer(pkfit_1cmp_tran3)
#Calculate BIC
bic(pkfit_1cmp_tran3)
 
####TESTING 1 CMP 4 TRAN FOR COMPLETENESS ##############
###ONE-COMPARTMENT 4 TRANS####
pk_1cmp_tran4 = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =4)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init = 45)
        """
        Mean transit time
        """
        tMTT ∈ RealDomain(; lower = 0, init = 5)
        """
          - ΩCL
          - ΩVc
          - ΩMTT
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [ 0.001, 0.001, 0.001])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.38)
    end

    @random begin
        η ~ MvNormal(Ω)
    end
    @pre begin
        CL = tvcl * exp(η[1])
        Vc = tvv * exp(η[2])
        MTT = tMTT * exp(η[3])

        ktr = 5/ MTT
        k20= CL/Vc
    end

    @dynamics begin
        Abs' = -ktr*Abs
        B1' = ktr *Abs- ktr*B1
        B2' = ktr* B1 - ktr*B2
        B3' = ktr* B2 -ktr*B3
        B4' = ktr* B3 -ktr*B4
        Central' = ktr*B4 - k20 * Central
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
pkparam1tran4 = init_params(pk_1cmp_tran4)
#Fit using FOCE
pkfit_1cmp_tran4 = fit(pk_1cmp_tran4, PYZ_pk, pkparam1tran4, FOCE())
#Inference
infer(pkfit_1cmp_tran4)
#Calculate BIC
bic(pkfit_1cmp_tran4)


#####################Testing a 2cmp model #####################################################
########################
#Testing 2 cmp 1 transit 
pk_2cmp_tran1 = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =5)
        """
        Clearance Periph (L/hr)
        """
        tvclp ∈ RealDomain(; lower = 0, init =3)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init = 40)
        """
        Volume Periph (L)
        """
        tvvp ∈ RealDomain(; lower = 0, init = 55)
        """
        Absorption rate
        """
        tMTT ∈ RealDomain(; lower = 0, init = 1)
        """
          - ΩCL
          - ΩVc
          - ΩCLp
          - ΩVp
          - Ωka
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.38)
    end

    @random begin
        η ~ MvNormal(Ω)
    end
    @pre begin
        CL = tvcl * exp(η[1])
        Vc = tvv * exp(η[2])
        CLp = tvclp * exp(η[3])
        Vp = tvvp * exp(η[4])
        MTT = (tMTT) * exp(η[5])

        ktr = 2/MTT
        k20= CL/Vc
        k23= CLp/Vc
        k32 = CLp/Vp
    end
    @dosecontrol begin
        bioav = ( ; Abs=1 * exp(η[6]))
    end
    @dynamics begin
        Abs' = -ktr*Abs
        B1'=ktr*Abs-ktr*B1
        Central' = ktr*B1 + k32 * Periph - k20 * Central
        Periph'= k23*Central - k32 *Periph
    end

    @derived begin
        cp := @. log(1+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p) 
    end

end
#Initial parms
pkparam_2tran = init_params(pk_2cmp_tran1)
#Fit using FOCE
pkfit_2cmp_tran1 = fit(pk_2cmp_tran1, PYZ_pk, pkparam_2tran, FOCE())
#Inference
infer(pkfit_2cmp_tran1)
#Calculate BIC
bic(pkfit_2cmp_tran1)

########################
#Testing 2 cmp 2 transit 
pk_2cmp_tran2 = @model begin
    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0, init =4)
        """
        Clearance Periph (L/hr)
        """
        tvclp ∈ RealDomain(; lower = 0, init =1)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0, init = 25)
        """
        Volume Periph (L)
        """
        tvvp ∈ RealDomain(; lower = 0, init = 50)
        """
        Absorption rate
        """
        tMTT ∈ RealDomain(; lower = 0, init = 1)
        """
          - ΩCL
          - ΩVc
          - ΩCLp
          - ΩVp
          - Ωka
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0, init = 0.38)
    end

    @random begin
        η ~ MvNormal(Ω)
    end
    @pre begin
        CL = tvcl * exp(η[1])
        Vc = tvv * exp(η[2])
        CLp = tvclp * exp(η[3])
        Vp = tvvp * exp(η[4])
        MTT = (tMTT) * exp(η[5])

        ktr = 3/MTT
        k20= CL/Vc
        k23= CLp/Vc
        k32 = CLp/Vp
    end
    @dosecontrol begin
        bioav = ( ; Abs=1 * exp(η[6]))
    end
    @dynamics begin
        Abs' = -ktr*Abs
        B1'=ktr*Abs-ktr*B1
        B2'=ktr*B1-ktr*B2
        Central' = ktr*B2 + k32 * Periph - k20 * Central
        Periph'= k23*Central - k32 *Periph
    end

    @derived begin
        cp := @. log(1+(Central / Vc)/1000) 
        """
        Concentration (ng/mL)
        """
        LOG_DV ~ @. Normal(cp, σ_p) 
    end

end
#Initial parms
pkparam_2tran2 = init_params(pk_2cmp_tran2)
#Fitting using FOCE
pkfit_2cmp_tran2 = fit(pk_2cmp_tran2, PYZ_pk, pkparam_2tran2, FOCE())
#Inference
infer(pkfit_2cmp_tran2)
#Calculating BIC
bic(pkfit_2cmp_tran2)
 