//
// Created by max on 2021-11-17.
//

#include "permuto_wrapper.h"

py::array_t<float> permuto_gauss_transform(py::array_t<float> sources, py::array_t<float>targets, py::array_t<float> weights, py::array_t<float> bandwidths)
{
    int d = sources.shape(1);
    int N = sources.shape(0);
    int R = targets.shape(0);

    auto rs = sources.unchecked<2>();
    auto rt = targets.unchecked<2>();
    auto rw = weights.unchecked<1>();
    auto rb = bandwidths.unchecked<1>();

    float * pSources = new float[N * d];
    for( int i=0; i<N; i++ ){
        for( int j=0; j<d; j++ ){
            pSources[d*i + j] = rs(i,j) / rb(j);
        }
    }

    float * pTargets = new float[R * d];
    for( int i=0; i<R; i++ ){
        for( int j=0; j<d; j++ ){
            pTargets[d*i + j] = rt(i,j) / rb(j);
        }
    }

    float * pWeights= new float[N];
    for( int i=0; i<N; i++ ){
        pWeights[i] = rw(i);
    }

    float * pOut = new float[R];

    Permutohedral p;
    p.init(pSources, pTargets, d, N, R);
    p.compute(pOut, pWeights, 1);

    py::array_t<float> output_values(R);
    auto ro = output_values.mutable_unchecked<1>();
    for( int i=0; i<R; i++ ){
        ro(i) = pOut[i];
    }

    delete[] pSources;
    delete[] pTargets;
    delete[] pWeights;
    delete[] pOut;

    return output_values;
}

