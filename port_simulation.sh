#!/usr/bin/bash

# Define the base command for running the ns-3 simulation
base_command="./ns3 run scratch/trail_error --"

# Define the CSV file for the simulations
csv_file="norm_off_port02502500.csv"

# Write the header if the file does not exist
if [ ! -f "$csv_file" ]; then
    echo "portPowerLabel,enableMimoFeedback,rankLimit,gnbUeDistance,Metric,Value" > "$csv_file"
fi

# Set fixed values as per requirements
enable_mimo_feedback=1  # MIMO feedback enabled by default
rank_limit=4            # Rank fixed to 4

# Loop through distances from 0 to 1000 with an increment of 20
for distance in $(seq 0 20 1000); do
    # Construct the full simulation command with default parameters
    full_command="$base_command --numColumnsGnb=2 --numRowsGnb=2 --numHPortsGnb=2 --numVPortsGnb=1 --xPolGnb=1 --numColumnsUe=2 --numRowsUe=1 --numHPortsUe=2 --numVPortsUe=1 --xPolUe=1 --fullSearchCb=ns3::NrCbTypeOneSp --enableMimoFeedback=$enable_mimo_feedback --rankLimit=$rank_limit --gnbUeDistance=$distance "
    
    # Notify user of current simulation
    echo "Running simulation with distance=$distance"
    
    # Run the simulation and capture output
    output=$(eval "$full_command")
    
    # Process simulation output and append to the CSV file
    echo "$output" | grep -E "Tx Packets:|Tx Bytes:|TxOffered: |Rx Bytes: |Throughput: |Mean delay: |Mean jitter: |Rx Packets: " | 
    awk -v lbl="default" -v emf="$enable_mimo_feedback" -v rl="$rank_limit" -v dist="$distance" '
    {
        if ($1 == "TxOffered:") {
            value = $2
            print lbl "," emf "," rl "," dist ",TxOffered," value
        }
        if ($1 == "Throughput:") {
            value = $2
            print lbl "," emf "," rl "," dist ",Throughput," value
        }
        if ($1 == "Tx" && $2 == "Packets:") {
            value = $3
            print lbl "," emf "," rl "," dist ",Tx Packets," value
        }
        if ($1 == "Tx" && $2 == "Bytes:") {
            value = $3
            print lbl "," emf "," rl "," dist ",Tx Bytes," value
        }
        if ($1 == "Rx" && $2 == "Bytes:") {
            value = $3
            print lbl "," emf "," rl "," dist ",Rx Bytes," value
        }
        if ($1 == "Rx" && $2 == "Packets:") {
            value = $3
            print lbl "," emf "," rl "," dist ",Rx Packets," value
        }
        if ($1 == "Mean" && $2 == "delay:") {
            value = $3
            print lbl "," emf "," rl "," dist ",Mean delay," value
        }
        if ($1 == "Mean" && $2 == "jitter:") {
            value = $3
            print lbl "," emf "," rl "," dist ",Mean jitter," value
        }
    }' >> "$csv_file"
done

echo "All simulations completed"
