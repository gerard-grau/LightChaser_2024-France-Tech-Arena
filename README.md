# Optical Network Service Continuity Optimization

This project contains our solution for the **Tech Arena 2024 - Optical Network Service Continuity Optimization Challenge**, developed by the team **Gerard Grau, Marc Herrero, and Pol Resina**.

## Problem Overview

The challenge simulates an **optical network** represented as an undirected connected graph with nodes and edges, where each edge corresponds to a real-world optical fiber containing multiple wavelengths. Multiple optical services are running through the network, each requiring specific wavelengths along their paths between a source and a destination.

The objective is to design an algorithm that **replans service paths** when fibers fail, ensuring maximum service continuity while adhering to several constraints related to wavelength allocation and channel conversion. Services that cannot be replanned are considered **dead**, and their resources remain occupied, making them unavailable for future planning.

### Key Constraints:
- Services must maintain their original bandwidth and path attributes.
- Wavelengths must be consistent across edges, except where channel conversion is available.
- A service can reuse its own resources from the old path during replanning but cannot use the resources from other services.
- Replanning is performed simultaneously for all affected services after a fiber failure.

The goal is to **maximize the total value of surviving services** after multiple fiber failures.

## Our Solution

Our solution focuses on ensuring service continuity while optimizing the allocation of wavelengths and making efficient use of channel conversion opportunities. Here's a breakdown of our approach:

### 1. **Graph Representation and Initialization**
We represent the optical network as an undirected graph, where each edge corresponds to an optical fiber with a defined set of wavelengths. We initialize the graph based on the input parameters, including the nodes, edges, services, and their respective wavelength allocations.

### 2. **Service Path Replanning Algorithm**
Upon the failure of a fiber (edge), services that traverse this edge need to be replanned. The main challenge is to find alternative paths for these services, ensuring that:
- The new paths have sufficient available wavelengths.
- Channel conversion opportunities are used efficiently, minimizing the number of conversions and ensuring wavelength continuity as much as possible.
  
We use a **modified Dijkstra’s algorithm** to find the shortest available path for each service. Our implementation optimizes for wavelength reuse along edges and ensures that channel conversion constraints are respected. The algorithm iterates through all affected services and tries to allocate the required contiguous wavelengths on each new path.

### 3. **Handling Resource Conflicts**
One of the most critical aspects of our solution is managing the allocation of wavelengths across different services. Services must not overlap in their wavelength usage, and services affected by the same fiber failure must compete for available resources. Our approach incorporates:
- **Wavelength Conflict Resolution**: Services are prioritized based on their value (i.e., the importance of keeping them alive), and we implement a conflict resolution mechanism that ensures high-value services get replanned first.
- **Efficient Channel Conversion**: Nodes with channel conversion capabilities are leveraged strategically to reassign wavelengths where continuity is broken, minimizing service downtime.

### 4. **Dynamic Replanning During Fiber Failures**
Whenever a fiber fault occurs, the services affected by the fault are dynamically replanned. This involves recalculating paths and checking wavelength availability for each service. The algorithm handles multiple services simultaneously, ensuring that replanning is done in parallel without interfering with other services' resources.

The replanning process involves:
- Identifying affected services after each failure.
- Recalculating the shortest path with available wavelengths for each service.
- Applying channel conversion where necessary to maintain wavelength continuity.

### 5. **Optimizing the Survival Value of Services**
The key metric for scoring our solution is the total value of surviving services after each test scenario. Our algorithm focuses on maximizing this value by:
- Prioritizing high-value services during the replanning phase.
- Efficiently using available wavelengths and channel conversion resources.
- Minimizing the number of services that remain "dead" after replanning.

### Performance Considerations
Our solution was designed to operate within the strict time (90 seconds) and memory (512MB) limits imposed by the challenge. We optimized the graph traversal and pathfinding algorithms to ensure scalability even with a large number of services and fiber failures.
