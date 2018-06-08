/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * Authors: "Jian Wang(jianopt@163.com)"
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#pragma once

#include "cppsci_macros.hpp"
#include "cppsci_array.hpp"
#include "cppsci_fft.hpp"
#include "cppsci_math.hpp"
#include "cppsci_signal.hpp"
#include "cppsci_time.hpp"

/////////////////////////////////////////////////////////////
// NOTE: These functions refer to MATLAB DeepLearnToolbox.
/////////////////////////////////////////////////////////////

namespace cppsci { namespace mlearn {

/**
 * Layer of convolution neural network.
 */
struct CNNLayer {
    char type;
    V<Arrayd> a;
    V<Arrayd> d;
    V<double> b, db;
    V<V<Arrayd>> k, dk;
    int scale;
    int outputmaps;
    int kernelsize;
};

/**
 * Construct a cnn input layer.
 */
CNNLayer cnn_input_layer();

/**
 * Construct a cnn convolution layer.
 */
CNNLayer cnn_conv_layer(int outputmaps_, int kernelsize_);

/**
 * Construct a cnn sampling layer.
 */
CNNLayer cnn_samp_layer(int scale_);

/**
 * CNN options.
 */
struct CNNOpts {
    int alpha;
    int batchsize;
    int numepochs;
};

/**
 * Colvolution neural network.
 */
struct CNN {
    V<CNNLayer> layers;
    Arrayd ffb, dffb, ffW, dffW, fv, fvd;
    Arrayd e, o, od;
    Vd rL;
    double L;
};

template<typename _V>
Arrayd sigm(const _V &v) {
    Arrayd m(v.shape());
    for (int i = 0; i < m.size(); i++) m[i] = 1.0 / (1.0 + STD_ exp(-v[i]));
    return STD_ move(m);
}

/**
 * Setup the cnn.
 */
void cnn_setup(CNN &net, const Arrayd &x, const Arrayd &y);

/**
 * CNN feedforward.
 */
void cnn_ff(CNN &net, Arrayd &x);

/**
 * CNN back propogation.
 */
Arrayd cnn_bp(CNN &net, Arrayd &y);

/**
 * CNN test.
 */
STD_ pair<double, Vi> cnn_test(CNN &net, Arrayd &x, Arrayd &y);

/**
 * CNN train.
 */
void cnn_train(CNN &net, Arrayd &x, Arrayd &y, const CNNOpts &opts);

/**
 * Apply gradients.
 */
void cnn_applygrads(CNN &net, const CNNOpts &opts);

/**
 * Save the CNN object to file.
 */
void cnn_save(const CNN &cnn, Str filename);

/**
 * Test the CNN using trained parameters.
 */
Arrayd cnn_test_m(CNN &cnn, Arrayd &x);

//void cnn_numgradcheck(CNN &net, Arrayd &x, Arrayd &y);

/**
 * Serialize a CNNLayer object to the stream.
 */
void serialize(STD_ ostream &stream, const CNNLayer &l);

/**
 * Parse a CNNLayer object from the stream.
 */
void parse(STD_ istream &stream, CNNLayer &l);

/**
 * Serialize a CNN object to the stream.
 */
void serialize(STD_ ostream &stream, const CNN &cnn);

/**
 * Parse a CNN object from the stream.
 */
void parse(STD_ istream &stream, CNN &cnn);

}} // namespace cppsci::mlearn

