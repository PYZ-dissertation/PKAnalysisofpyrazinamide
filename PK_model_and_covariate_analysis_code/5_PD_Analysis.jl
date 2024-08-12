############
#Julia code to perform some tests on Bacterial clearance including logistic regression models
# and Weibull time to event models. Some plots are produced to compare basic differences in bacterial clearance
#between HIV-negative and HIV positive patients
############

#Load packages
using Pkg, XLSX, Pumas, Dates, DataFrames, DataFramesMeta, DataStructures, GLM, Survival, Plots, PumasUtilities, Pumas.Latexify,AlgebraOfGraphics,
Makie, StatsPlots
using AlgebraOfGraphics: data as AoG_data, mapping, visual, draw, Axis

#Set working directory
######
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
    if col in ["AMT", "AGE", "ALTRES", "TOTBRES", "CREATRES", "WT"]
        PYZ_data[!, col] = convert(Vector{Float64}, PYZ_data[!, col])
    elseif col in ["CMT", "EVID", "ID", "MDV", "SEX", "HIV"]
        PYZ_data[!, col] = convert(Vector{Int64}, PYZ_data[!, col])
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

# Apply the function to each group and update the TIME column
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

#Sort data by IDs then TAD
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

#Add column for time as a covariate
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

#Function to calculate AUC for each individual
function calculate_auc(time_points, dv)
    auc = 0.0
    for i in 1:(length(time_points) - 1)
        dt = time_points[i+1] - time_points[i]
        auc += ((dv[i]) + (dv[i+1])) / (2*1000) * dt # mg/L
    end
    return auc
end

# Function to calculate CMAX
function calculate_cmax(dv)
    return maximum(dv)
end

# Initialize a dictionary to store AUC results and CMAX for each ID
auc_results = OrderedDict{Int, Vector{Float64}}()
cmax_results = OrderedDict{Int, Float64}()

# Group data by subjects and TP2 and calculate AUC for each combination
cleaned_data = dropmissing(PYZ_data)
grouped_data = groupby(cleaned_data, [:ID, :TP2])

for (key, group) in pairs(grouped_data)
    id, tp2_status = key
    if tp2_status in [0, 1]
        # Calculate AUC for the subject
        auc = calculate_auc(group.TAD, group.DV)
        
        # Store the AUC in the corresponding ID entry
        if haskey(auc_results, id)
            push!(auc_results[id], auc)
        else
            auc_results[id] = [auc]
        end
        
        # Calculate CMAX for the subject and store it (CMAX in whole trial)
        if !haskey(cmax_results, id)
            cmax_results[id] = (calculate_cmax(group.DV)/1000)
        end
    end
end

# Create a new DataFrame with the average AUC and CMAX for each ID
ids = collect(keys(auc_results))
avg_auc = [mean(auc_results[id]) for id in ids]
cmax = [cmax_results[id] for id in ids]

results_df = DataFrame(ID = ids, AUC = avg_auc, CMAX = cmax)

#Loading Pyrazinamide data
PD_data = DataFrame(XLSX.readtable("PD_data_full.xlsx", 1, infer_eltypes = true))

#Change time from weeks to hours
PD_data.TIME = (PD_data.Bact_clearance)*24*7

#Replace missing time (those who didnt complete trial) with 0
replace!(PD_data.TIME, missing => 0);

#Create DV column, 0 if Bact_clearance = missing (i.e. never cleared), 1 otherwise
PD_data.DV = ifelse.(ismissing.(PD_data.Bact_clearance), 0, 1)

#EVID column
PD_data.EVID = ifelse.(isinteger.(PD_data.ID), 0, 1)

#Add AMT column for PD functions
@rtransform! PD_data :AMT = 0; 

#Add AUC and CMAX data too PD_data
PD_data = innerjoin(PD_data, results_df, on = :ID)

#Change miss spelling of Thailand
replace!(PYZ_data.District, "THAILNAD" => "THAILAND")

#Access district of each ID
unique_ids_df = unique(PYZ_data[:, [:ID, :District, :SEX]])

# Join the unique IDs DataFrame with the new DataFrame to add District information
PD_data = leftjoin(PD_data, unique_ids_df[:, [:ID, :District,:SEX]], on = :ID)

#Creating a flag column if below reccomended minimum AUC of 363
PD_data.AUC_flag = ifelse.(PD_data.AUC .< 363, 0, 1)

#data now in correct format
##############

#Population with HIV as covariate- HIV analysis
PD_HIV = read_pumas(
    PD_data;
    observations = [:DV],
    id = :ID,
    time = :TIME,
    evid = :EVID,
    amt = :AMT,
    covariates = [:HIV]
    )

#Population with SEX as a covariate
PD_SEX = read_pumas(
    PD_data;
    observations = [:DV],
    id = :ID,
    time = :TIME,
    evid = :EVID,
    amt = :AMT,
    covariates = [:SEX]
    )

#Dropping those who were not cured 
Delay_PD = @chain PD_data begin
    dropmissing(:DelayBC)
end

 
#Logistic Regressions
# Fit the logistic regression model
#HIV logistic regression model
formula_HIV = @formula(DelayBC ~ HIV )
logit_model_HIV = glm(formula_HIV, PD_data, Binomial(), LogitLink())

#Age logistic regression model
formula_AGE= @formula(DelayBC ~ AGE )
logit_model_AGE = glm(formula_AGE, PD_data, Binomial(), LogitLink())

#AUC logistic regression model
formula_AUC = @formula(DelayBC ~ AUC)
logit_model_AUC = glm(formula_AUC, PD_data, Binomial(), LogitLink())

#CMAX logistic regression model
formula_CMAX = @formula(DelayBC ~ CMAX)
logit_model_CMAX = glm(formula_CMAX, PD_data, Binomial(), LogitLink())

#WTM0 logistic regression model
formula_WtM0= @formula(DelayBC ~ WtM0)
logit_model_WtM0 = glm(formula_WtM0, PD_data, Binomial(), LogitLink())

#SEX logistic regression model
formula_SEX = @formula(DelayBC ~ SEX )
logit_model_SEX = glm(formula_SEX, PD_data, Binomial(), LogitLink())

#AUC_flag logistic regression model
formula_AUC_flag = @formula(DelayBC ~ AUC_flag )
logit_model_AUC_flag= glm(formula_AUC_flag, PD_data, Binomial(), LogitLink())

#Time to event Weiburg model with HIV as a covariate
PD_HIVmod = @model begin
        @param begin
            λ₁ ∈ RealDomain(; lower = 0, init = 0.001)
            β ∈ RealDomain(; init = 0.001)
            Κ ∈ RealDomain(; lower = 0, init = 0.001)
        end
    
        @covariates HIV
    
        @pre begin
            _Κ = Κ
            _λ₀ = λ₁ * exp(β * HIV) #λ₀ = total hazard
        end
    
        @vars begin
            # Weibull
            # 1e-10 for model numerical stability
            λ = _λ₀ * _Κ * (_λ₀ * t + 1e-10)^(_Κ - 1) #λ = PDF of the Weibull distribution 
        end
    
        @dynamics begin
            # the derivative of cumulative hazard is equal to hazard
            Λ' = λ #Λ  = cumulative hazard
        end
    
        @derived begin
            # Pumas derived TTE function
            DV ~ @. TimeToEvent(λ, Λ)
        end
end

#Time to event Weiburg model with SEX as a covariate
PD_SEXmod = @model begin
    @param begin
        λ₁ ∈ RealDomain(; lower = 0, init = 0.001)
        β ∈ RealDomain(; init = 0.001)
        Κ ∈ RealDomain(; lower = 0, init = 0.001)
    end

    @covariates SEX

    @pre begin
        _Κ = Κ
        _λ₀ = λ₁ * exp(β * SEX) #λ₀ = total hazard
    end

    @vars begin
        # Weibull
        # 1e-10 for model numerical stability
        λ = _λ₀ * _Κ * (_λ₀ * t + 1e-10)^(_Κ - 1) #λ = PDF of the Weibull distribution 
    end

    @dynamics begin
        # the derivative of cumulative hazard is equal to hazard
        Λ' = λ #Λ  = cumulative hazard
    end

    @derived begin
        # Pumas derived TTE function
        DV ~ @. TimeToEvent(λ, Λ)
    end
end
#fit time to event Wiebull model for HIV covariate 
PD_HIV_fit = fit(PD_HIVmod, PD_HIV, init_params(PD_HIVmod), Pumas.NaivePooled());
#fit time to event wiebull model for SEX covariate 
PD_SEX_fit = fit(PD_SEXmod, PD_SEX, init_params(PD_SEXmod), Pumas.NaivePooled());

#Produce vpc plots
vpc_pd_final = vpc(
        PD_HIV_fit; 
        samples=2000,
        covariates = [:time],
        stratify_by = [:HIV],
        smooth= true
        )

vpc_plot(
    PD_HIVmod,
    vpc_pd_final;
    figure = (; resolution = (1400, 1000), fontsize = 22),
    axis = (;
        xlabel = "Time after dose (hr)",
        ylabel = "Observed/Predicted\n Log Concentration (ng/mL)",
    ),
    facet = (; combinelabels = true),
)

#Dataframes of fits
PD_wei_HIV_df = (infer(PD_HIV_fit)) #parameter fits for HIV model
PD_wei_SEX_df = (infer(PD_SEX_fit)) #parameter fits for SEX model

#The below functions use the parameter output values from the fit objects
#They generate the fitted curve (as is generated by the fit objects)
#create function for post-processing HIV
function hazard_weibull_HIV(λ₁, β, Κ, HIV, T)
    λ₀ = λ₁ * exp(β * HIV)                           # estimated hazard
    λ = [λ₀ * Κ * (λ₀ * t + 1e-10)^(Κ - 1) for t = 1:T] # cumulative hazard
    return 1 .- cumsum(λ)
end

#create function for post-processing SEX
function hazard_weibull_SEX(λ₁, β, Κ, SEX, T)
    λ₀ = λ₁ * exp(β * SEX)                          # estimated hazard
    λ = [λ₀ * Κ * (λ₀ * t + 1e-10)^(Κ - 1) for t = 1:T] # cumulative hazard
    return 1 .- cumsum(λ)
end

#Time seriesof 14 weeks
t = 1:(168*10)

#Below, the hazard_weibull function is called to generate the fitted data points
#HIV
λ_weibull_HIV_0 = hazard_weibull_HIV(
    0.00074237, # estimated λ₁
    0.23383,  # estimated β
    1.7063,    # Κ 
    0,          # HIV value
    length(t),  # T
)

λ_weibull_HIV_1 = hazard_weibull_HIV(
    0.00074237, # estimated λ₁
    0.23383,  # estimated β
    1.7063,    # Κ 
    1,          # HIV value
    length(t),  # T
)

#SEX
λ_weibull_SEX_0 = hazard_weibull_SEX(
    0.00074259, # estimated λ₁
    0.22662  ,  # estimated β
    1.7282,    # Κ 
    0,          # SEX value
    length(t),  # T
)

λ_weibull_SEX_1 = hazard_weibull_SEX(
    0.00074259, # estimated λ₁
    0.22662  ,  # estimated β
    1.7282,    # Κ 
    0,          # SEX value
    length(t),  # T
)

#Plot the HIV covariate model
plt_hiv =
    AoG_data(
    # data in a long table format
        (;
        t = repeat(t, 2),             # 2 DOSES ⋅ 1 AFT models
        λs = vcat(
            λ_weibull_HIV_0,
            λ_weibull_HIV_1,
        ),
        HIV = repeat(
            [0, 1];
            inner = length(t), # repeat 0 500x then 1 500x
            outer = 1,          # 1 AFT models
        ),
        MODEL = repeat(["Weibull"]; inner = 2 * length(t)),    # repeat each model 1_000x
    )) *
    mapping(:t => L"t", :λs => L"S(t)"; color = :HIV => nonnumeric) *
    visual(Lines)
draw(plt_hiv; axis = (; xticks = 0:168:1680, yticks = 0.0:0.25:1.0, limits = ((0, nothing), (0, nothing))))

#Plot the SEX covariate model
plt_SEX =
    AoG_data(
    # data in a long table format
        (;
        t = repeat(t, 2),             # 2 DOSES ⋅ 1 AFT models
        λs = vcat(
            λ_weibull_SEX_0,
            λ_weibull_SEX_1,
        ),
        SEX = repeat(
            [0, 1];
            inner = length(t), # repeat 0 500x then 1 500x
            outer = 1,          # 1 AFT models
        ),
        MODEL = repeat(["Weibull"]; inner = 2 * length(t)),    # repeat each model 1_000x
    )) *
    mapping(:t => L"t", :λs => L"S(t)"; color = :SEX => nonnumeric) *
    visual(Lines)
draw(plt_SEX; axis = (; xticks = 0:168:1680, yticks = 0.0:0.25:1.0, limits = ((0, nothing), (0, nothing))))


#Kaplan mieir step function
tps = PD_data.TIME #time series
clr = PD_data.DV #event occurrance

survival_fit = fit(KaplanMeier, tps, clr) #fit plot

#draw plot
plt =
    AoG_data((; t = survival_fit.events.time, s_t = survival_fit.survival)) *
    mapping(:t => L"t", :s_t => L"S(t)") *
    visual(Stairs; step = :post, linewidth = 3)
draw(
    plt;
    axis = (; xticks = 1:336:3700, yticks = 0.0:0.1:1.0, limits = ((0, nothing), (0, nothing))),
)

# Create data for the plot
plot_data = DataFrame(t = survival_fit.events.time, s_t = survival_fit.survival)

#Making some basic diagnostic plots#
# Sort DataFrame by enrolment_date
sort!(PD_data, :enrolment_date)

#Change HIV status to binary
PD_data.hivpost = ifelse.(PD_data.HIV .== 0, "Negative", "Positive")
PD_data2 = dropmissing(PD_data)

# Create grouped bar chart
@df PD_data2 groupedbar( :Bact_clearance, group = :hivpost,
            xlabel = "Bacterial Clearance Time (Weeks)", ylabel = "Frequency",
            title = "Grouped Bar Chart of Bacterial Clearance Time by HIV Status")

#Create a grouped histogram           
@df PD_data2 groupedhist( :Bact_clearance, group = :hivpost,
            xlabel = "Bacterial Clearance Time (Weeks)", ylabel = "Frequency",
            title = "Grouped Histogram of Bacterial Clearance Time by HIV Status")

# Box plot 
Plots.plot(
    Plots.boxplot([PD_data2[PD_data2.hivpost .== "Negative", :Bact_clearance], PD_data2[PD_data2.hivpost .== "Positive", :Bact_clearance]],
            boxalpha = 0.8,  # Adjust box transparency
            boxcolor = [:blue, :green],  # Specify box colors
            whisker_color = :black,  # Specify whisker color
            marker_size = 4,  # Adjust marker size
            legend = false),  # No legend
    xlabel = "HIV Status",  # X-axis label
    ylabel = "Bacterial Clearance Time (Weeks)",  # Y-axis label
    title = "Box Plot of Bacterial Clearance Time by HIV Status",  # Plot title
    xticks = (1:2, ["Negative", "Positive"]),  # X-axis ticks and labels
    linecolor = :black,  # Specify line color
    linewidth = 1.5,  # Adjust line width
    framestyle = :box,  # Add a box around the plot
    grid = false  # Disable grid
)

#Print some summary stats
@by(PD_data2, :hivpost,
    :mean_Bact_clearance = mean(:Bact_clearance),
    median_Bact_clearance = median(:Bact_clearance),
    std_Bact_clearance = std(:Bact_clearance),
    min_Bact_clearance = minimum(:Bact_clearance),
    max_Bact_clearance = maximum(:Bact_clearance),
    q25_Bact_clearance = quantile(:Bact_clearance, 0.25),
    q75_Bact_clearance = quantile(:Bact_clearance, 0.75)
)
