# NS3 Port Power Allocation Simulation

This guide walks you through setting up and running the NS3 port power allocation simulation.

## Prerequisites
- NS3 development environment properly installed
- Access to the `port_power_allocation.cc` file which is in the github repo 
- Basic familiarity with terminal/command line
- For detailed NS3 NR module documentation, visit: [NS3 NR Documentation](https://cttc-lena.gitlab.io/nr/html/)
- we recommend using ubuntu 20.04 and 5G Lena module 3.3 and Ns3 version 3.42

## Setup and Execution

### Step 1: File Placement
Navigate to the NS3 development directory and place the scratch folder :
```bash
cd ns3-dev
cd scratch
```

### Step 2: Configure Port Power Parameters
Edit the `port_power_allocation.cc` file to configure port parameters. **Important**: Ensure the number of ports in your power allocation parameter matches the total ports in your scenario.

**Port Calculation Formula:**
- If `dual_polarization = true`: Total ports = (num_horizontal_ports × num_vertical_ports) × 2
- If `dual_polarization = false`: Total ports = num_horizontal_ports × num_vertical_ports

**Example:**
```cpp
num_horizontal_ports = 1
num_vertical_ports = 2  
dual_polarization = true
// Total ports = (1 × 2) × 2 = 4 ports
```

Make sure your port power allocation parameter has exactly the same number of ports as the simulation parameters.

### Step 3: Save and Build
Save your changes and build the simulation:
```bash
cd ns3-dev
./ns3 build
```

### Step 4: Run Single Scenario
Execute the simulation:
```bash
./ns3 run scratch/port_power_allocation
```

Results will be displayed in the terminal.

## Running Multiple Scenarios (Graph Generation)

### Step 5: Setup Multiple Scenarios
Repeat steps 1-3 for each scenario you want to include in your graph.

### Step 6: Execute Batch Simulation
Make the simulation script executable and run it:
```bash
chmod +x port_simulation.sh
./port_simulation.sh
```

**Note:** Make sure you run these commands from within the `ns3-dev` directory.

## Output
- Single scenario results: Displayed in terminal
- Multiple scenario results: Saved as CSV files in the `ns3-dev` directory
- Generated graphs: Also saved in the `ns3-dev` directory

## Troubleshooting
- Ensure all commands are executed from the `ns3-dev` directory
- Verify that the number of ports in your power allocation matches your scenario configuration
- Check that all required files are present before building
- Make sure NS3 is properly installed and configured

## File Structure
```
ns3-dev/
├── scratch/
│   └── port_power_allocation.cc
├── port_simulation.sh
└── [output files will be generated here]
```

## Contact
For any inquiries, please don't hesitate to reach out:
```
tarekomar492@gmail.com
```

## Additional Resources
- [NS3 NR Module Documentation](https://cttc-lena.gitlab.io/nr/html/) - Comprehensive documentation for the NS3 NR module

---
**Thanks**
