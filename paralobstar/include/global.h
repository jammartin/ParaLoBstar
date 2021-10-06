//
// Created by Johannes Martin on 20.09.21.
//

#ifndef PARALOBSTAR_GLOBAL_H
#define PARALOBSTAR_GLOBAL_H

namespace global {
    constexpr int dim { 3 };
    constexpr int powdim { 1 << dim }; // 2^dim
    constexpr int maxTreeLvl { 64/dim }; // using 64-bit key
    constexpr double G { 1. }; // gravitational constant
}

#endif //PARALOBSTAR_GLOBAL_H
