# Sketch-Arch
Sketch-Arch, an energy-efficient architecture for real-time mapping of long read sequences in portable systems. We present an implementation that targets the proposed PNM architecture, and also a GPU-friendly implementation for baseline comparisons. The implementations use a structured hash-map and can leverage at least two different types of 𝑘-mer sampling methods (minimizers and syncmers) under the hood. 

<p align="center">
  <img src="Sketch-Arch.png" alt="Image description" width="400"/>

</p>

Experimental results show that Sketch-Arch reduces memory footprint by up to 7.4× through sketch-based data reduction.

## Dependencies

gcc 8.5

nvcc 11.7

## Build
1. **Run:**

   ```bash
   make

## Execute
1. **Run:**

   ```bash
   ./script.sh
   

