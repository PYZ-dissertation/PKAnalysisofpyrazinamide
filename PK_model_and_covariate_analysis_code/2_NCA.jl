############
#Julia code to conduct two NCAs for the concentration data on day zero and week 6
############

#Load packages
using XLSX, Pumas, NCA, Dates, DataFrames, DataFramesMeta, NCAUtilities

######
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
# Add TP column to PYZ_data
PYZ_data = add_tp_column!(PYZ_data)

# Set ROUTE column to "ev" for all rows
PYZ_data.ROUTE .= "ev"

# Filter data for D0 and W6 time points
PYZ_dataD0 = filter(row -> occursin("D0", row."Time.point"), PYZ_data)
PYZ_dataW6 = filter(row -> occursin("W6", row."Time.point"), PYZ_data)

# Sort data by ID and TAD
PYZ_dataD0 = sort(PYZ_dataD0, [:ID, :TAD])
PYZ_dataW6 = sort(PYZ_dataW6, [:ID, :TAD])

# Multiply DV by 1000 for D0 and W6 data
PYZ_dataD0.DV .= PYZ_dataD0.DV*1000 
PYZ_dataW6.DV .= PYZ_dataW6.DV*1000 

#data now in correct format
##############

####NON-COMPARTMENTAL ANALYSIS#####
#NCA OF D0###########################################################################################################
#Creating an  NCA Populationfor day 0
PYZ_ncaD0 = read_nca(
    PYZ_dataD0;
    id = :ID,
    time =:TAD,
    amt = :AMT,
    observations = :DV,
    route = :ROUTE,
)

#Constructing a full NCA report for D0
PYZ_nca_repD0 = run_nca(PYZ_ncaD0;
                        studytitle ="PYZ D0",
                        author = [("Emily Hatchwell", "Pumas-AI")],
                        sigdigits = 3)
report(PYZ_nca_repD0 )

############NCA OF W6##############################################################################
#Creating an  NCA Population for week 6
PYZ_ncaW6 = read_nca(
    PYZ_dataW6;
    id = :ID,
    time =:TAD,
    amt = :AMT,
    observations = :DV,
    route= :ROUTE
)

#Constructing a full NCA Report for W6
PYZ_nca_repW6 = run_nca(PYZ_ncaW6;
                        studytitle ="PYZ W6",
                        author = [("Emily Hatchwell", "Pumas-AI")],
                        sigdigits = 3)
report(PYZ_nca_repW6 )

