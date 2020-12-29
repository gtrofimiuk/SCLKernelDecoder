# SCLKernelDecoder
This is an implementation of the window processing algorithm for `K_16` and `K_32` kernels from https://arxiv.org/abs/2010.07884. Kernel processors are embedded in [SCL](https://ieeexplore.ieee.org/document/7055304) decoder. This implementation also supports the decoding of [polar subcodes](https://ieeexplore.ieee.org/document/7339451). Simulations run in AWGN channel

## Build

The sources can be built with Visual Studio 2017

## Running
### Usage

The simulation program runs from the command line:
```
SCLKernelDecoder code_specification_file ListSize Eb/N0 MaxErrors MaxIterations
```
* `MaxIterations` is a maximal number of iterations performed by the decoder
* The simulation continues until `MaxErrors` decoding errors occurs

### Example usage
```
SCLKernelDecoder.exe 1024_512_Trofimiuk32_342.mpec 8 1.75 100 1000000
```

### Output line
`
Eb/N0 noise_standard_deviation block_error_rate ML_block_error_rate
`

### Code specification format
Here we introduce the code specification file format for polar subcodes with mixed kernels
```
<Length> <Dimension> <Minimum distance> <Number of layers> <Number of shortened symbols> <Number of shortened symbols> // Minimum distance can be set to zero
<Kernel_Name_1> <Kernel_Name_2> ... <Kernel_Name_lambda> //lambda is a number of layers

// R = Length - Dimension dynamic frozen symbols (DFS)
// The DFS is the constraint of the form u_i = u_{j_1} + u_{j_2} + u_{j_3} + ... + u_{j_w}, where all j_r < i
<w+1> <j_1> <j_2> ... <j_w> <i>
// Static frozen symbol i are represented as 
// 1 i
```

#### Example code specification
```
1024 1020 0 2 0 0
Trofimiuk32_342 Trofimiuk32_342

1 0
1 1
1 2
3 3 5 9
```

#### Supported kernels are 
* `Trofimiuk16_345` - 16 x 16 kernel `K_16`
* `Trofimiuk32_342` - 32 x 32 kernel `K_32`

This implementation does not support mixed kernels, shortening and puncturing. We introduce these parameters for futher implementations
