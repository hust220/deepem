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

#include "cppsci_mlearn.hpp"
#include "cppsci_serial.hpp"

namespace cppsci { namespace mlearn {

template<typename T>
static void print_array(T &&t) {
    int rows = t.shape(0);
    int cols = t.shape(1);
    std::cout << rows << ' ' << cols << std::endl;
    auto it = t.data();
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            std::cout << *it << ' ';
            it++;
        }
        std::cout << std::endl;
    }
}

CNNLayer cnn_input_layer() {
    CNNLayer l;
    l.type = 'i';
    return STD_ move (l);
}

CNNLayer cnn_conv_layer(int outputmaps_, int kernelsize_) {
    CNNLayer l;
    l.type = 'c';
    l.kernelsize = kernelsize_;
    l.outputmaps = outputmaps_;
    return STD_ move(l);
}

CNNLayer cnn_samp_layer(int scale_) {
    CNNLayer l;
    l.type = 's';
    l.scale = scale_;
    return STD_ move(l);
}

static Arrayd randmat(const Shape &shape) {
    Arrayd r(shape);
    Rand rand;
    for (int i = 0; i < r.size(); i++) {
        r[i] = rand((i+1)*100);
    }
    return std::move(r);
}

void cnn_setup(CNN &net, const Arrayd &x, const Arrayd &y) {
    int inputmaps = 1;
    Shape mapsize = {x.shape(1), x.shape(2)};
    JN_INFOV(mapsize);

    for (auto && l : net.layers) {
        JN_INFO(l.type);
        if (l.type == 's') {
            JN_INFO(inputmaps);
            for (auto && i : mapsize) i /= l.scale;
            l.b.resize(inputmaps, 0);
        }
        else if (l.type == 'c') {
            for (auto && i : mapsize) i -= l.kernelsize - 1;
            int fan_in = inputmaps * (l.kernelsize) * (l.kernelsize);
            int fan_out = (l.outputmaps) * (l.kernelsize) * (l.kernelsize);

            JN_INFO(inputmaps);
            JN_INFO(l.outputmaps);
            JN_INFO(l.kernelsize);

            l.b.resize(l.outputmaps, 0);

            l.k.resize(inputmaps);
            for (auto && kk : l.k) {
                kk.resize(l.outputmaps);
                for (auto && ll : kk) ll = Arrayd({l.kernelsize, l.kernelsize});
            }

#ifndef RAND
            Rand rand;
#endif
            for (int j = 0; j < l.outputmaps; j++) {
                for (int i = 0; i < inputmaps; i++) {
                    double *it = l.k[i][j].data();
                    for (int m = 0; m < l.kernelsize; m++) {
                        for (int n = 0; n < l.kernelsize; n++) {
#ifdef RAND
                            *it = (rand() - 0.5) * 2 * STD_ sqrt(6.0 / (fan_in + fan_out));
#else
                            *it = (rand((m*l.kernelsize+n+1)*100) - 0.5) * 2 * STD_ sqrt(6.0 / (fan_in + fan_out));
#endif
                            it++;
                        }
                    }
                    JN_INFO(i);
                    JN_INFO(j);
                    JN_INFO(fan_in);
                    JN_INFO(fan_out);
                    JN_INFOA(l.k[i][j]);
                }
            }
            inputmaps = l.outputmaps;
            JN_INFO(inputmaps);
        }
        else if (l.type == 'i') {
            // pass
        }
        else {
            // pass
        }
    }

    JN_INFOV(mapsize);
    JN_INFO(inputmaps);
    int fvnum = STD_ accumulate(mapsize.begin(), mapsize.end(), 1, [](int a, int b){return a*b;}) * inputmaps;
    JN_INFO(fvnum);
    int onum = y.shape(1);
    JN_INFO(onum);
    net.ffb = Arrayd({1, onum}, 0);
#ifdef RAND
    net.ffW = STD_ move(rand({fvnum, onum}).selfMinus(0.5).selfTimes(2.0 * STD_ sqrt(6.0 / (onum + fvnum))));
#else
    net.ffW = STD_ move(randmat({fvnum, onum}).selfMinus(0.5).selfTimes(2.0 * STD_ sqrt(6.0 / (onum + fvnum))));
#endif
    JN_INFOA(net.ffW);
}

void cnn_ff(CNN &net, Arrayd &x) {
    JN_INFOA(x);

    int n = net.layers.size();
    net.layers[0].a.resize(1);
    net.layers[0].a[0] = STD_ move(x);

    JN_INFOA(net.layers[0].a[0]);

    int inputmaps = 1;
    for (int l = 1; l < n; l++) {
        auto &layer = net.layers[l];
        auto &prev_layer = net.layers[l-1];
        if (layer.type == 'c') {
            JN_INFO(inputmaps);
            JN_INFO(layer.outputmaps);
            // below can probably be handled by insane matrix operations
            layer.a.resize(layer.outputmaps);

            Shape shapez = prev_layer.a[0].shape();
            JN_INFO(shapez[1]);
            JN_INFO(shapez[2]);
            shapez[1] -= layer.kernelsize-1;
            shapez[2] -= layer.kernelsize-1;
            JN_INFO(layer.kernelsize);
            JN_INFO(shapez[1]);
            JN_INFO(shapez[2]);

//#ifdef USEMPI
//            auto mpialloc = mpi_alloc(layer.outputmaps);
//#pragma omp parallel for
//            for (int j = mpialloc.beg; j < mpialloc.end; j++) {
//#else
#pragma omp parallel for
            for (int j = 0; j < layer.outputmaps; j++) {
//#endif
                Arrayd z(shapez, 0);
                for (int i = 0; i < inputmaps; i++) {
                    JN_INFO(i);
                    JN_INFOA(prev_layer.a[i]);
                    JN_INFOS("here3");
                    JN_INFOA(layer.k[i][j]);
                    z.selfPlus(signal::convolve(prev_layer.a[i], layer.k[i][j], signal::CONVOLVE_VALID, signal::CONVOLVE_PLAIN));
                }
                //  add bias, pass through nonlinearity
                layer.a[j] = Arrayd(shapez);
                double *a_it = layer.a[j].data();
                double *z_it = z.data();
                double b = layer.b[j];
                for (int i = 0; i < layer.a[j].size(); i++) {
                    *a_it = 1.0 / (1.0 + STD_ exp(-(*z_it + b)));
                    a_it++;
                    z_it++;
                }
                JN_INFOS("here1");
                JN_INFOA(layer.a[j]);
                //layer.a[j] = sigm(z + layer.b[j]);
            }
//#ifdef USEMPI
//            for (int j = 0; j < mpialloc.size; j++) {
//                if (mpi_rank() != j/mpialloc.bin) layer.a[j] = Arrayd(shapez);
//                JN_INFOS(string::format("Bcast from node %d with size of %d", j/mpialloc.bin, layer.a[j].size()));
//                JN_INFOS(string::format("total: %d, j: %d, Bin: %d, Node: %d, Size: %d", mpialloc.size, j, mpialloc.bin, mpi_rank(), layer.a[j].size()));
//                MPI_Bcast(layer.a[j].data(), layer.a[j].size(), MPI_DOUBLE, j/mpialloc.bin, MPI_COMM_WORLD);
//                JN_INFOS("Bcast done");
//            }
//#endif
            //  set number of input maps to this layers number of outputmaps
            inputmaps = layer.outputmaps;
        }
        else if (layer.type == 's') {
            JN_INFOS("downsample");
            JN_INFO(inputmaps);
            layer.a.resize(inputmaps);
            JN_INFO(layer.a.size());
            // downsample
            Arrayd ascale({layer.scale, layer.scale}, 1.0 / square(layer.scale));
            int scl = layer.scale;

//#ifdef USEMPI
//            Shape shape = signal::convolve(prev_layer.a[0], ascale, signal::CONVOLVE_VALID, signal::CONVOLVE_PLAIN).shape();
//            shape[1] = shape[1]/scl+1;
//            shape[2] = shape[2]/scl+1;
//            auto mpialloc = mpi_alloc(inputmaps);
//#pragma omp parallel for
//            for (int j = mpialloc.beg; j < mpialloc.end; j++) {
//#else
#pragma omp parallel for
            for (int j = 0; j < inputmaps; j++) {
//#endif
                JN_INFO(j);
                auto && z = signal::convolve(prev_layer.a[j], ascale, signal::CONVOLVE_VALID, signal::CONVOLVE_PLAIN);
                JN_INFOA(z);

                int zslis = z.shape(0);
                int zrows = z.shape(1);
                int zcols = z.shape(2);
                int zarea = zrows * zcols;

                int slis = zslis;
                int rows = zrows/scl + 1;
                int cols = zcols/scl + 1;

                layer.a[j] = Arrayd({slis, rows, cols});
                double *aj_it = layer.a[j].data();
                double *z_it = z.data();
                for (int sli = 0; sli < slis; sli++) {
                    double *z_row_it = z_it;
                    for (int row = 0; row < rows; row++) {
                        double *z_col_it = z_row_it;
                        for (int col = 0; col < cols; col++) {
                            *aj_it = *z_col_it;
                            aj_it++;
                            z_col_it += scl;
                        }
                        z_row_it += scl*zcols;
                    }
                    z_it += zarea;
                }

                JN_INFOS("here2");
                JN_INFOA(layer.a[j]);
            }
//#ifdef USEMPI
//            for (int j = 0; j < mpialloc.size; j++) {
//                if (layer.a[j].size() == 0) layer.a[j] = Arrayd(shape);
//                JN_INFOS(string::format("Bcast from node %d with size of %d", j/mpialloc.bin, layer.a[j].size()));
//                MPI_Bcast(layer.a[j].data(), layer.a[j].size(), MPI_DOUBLE, j/mpialloc.bin, MPI_COMM_WORLD);
//                JN_INFOS("Bcast done");
//            }
//#endif
        }
    }

    // concatenate all end layer feature maps into vector
    JN_INFOS("concatenate all end layer feature maps into vector");
    auto &layer = net.layers.back();
    JN_INFO(layer.a.size());

    int rows = layer.a[0].shape(0);
    int cols = 0;
    for (int j = 0; j < layer.a.size(); j++) {
        cols += layer.a[j].shape(1)*layer.a[j].shape(2);
    }
    JN_INFO(rows);
    JN_INFO(cols);

    net.fv = Arrayd({rows, cols});
    double *fv_it = net.fv.data();
    for (int j = 0; j < layer.a.size(); j++) {
        auto && sa = layer.a[j].shape();
//        JN_INFOV(sa);
//        JN_INFOA(layer.a[j]);
        int sz = sa[1]*sa[2];
        double *fv_it2 = fv_it;
        double *aj_it = layer.a[j].data();
        for (int k = 0; k < sa[0]; k++) {
            double *fv_it3 = fv_it2;
            for (int l = 0; l < sz; l++) {
                *fv_it3 = *aj_it;
                fv_it3++;
                aj_it++;
            }
            fv_it2 += cols;
        }
        //net.fv = hstack<double>(net.fv, reshape(layer.a[j], {sa[0], sa[1]*sa[2]}));
        fv_it += sz;
    }

    // feedforward into output perceptrons
    JN_INFOS("feedforward into output perceptrons");
    JN_INFOA(net.ffW);
    JN_INFOA(net.fv);
    JN_INFOA(net.ffb);
    JN_INFOA(multiply<double>(net.fv, net.ffW));
    //print_array(net.fv);
    //print_array(net.ffW);
    net.o = sigm(multiply<double>(net.fv, net.ffW) + repmat(net.ffb, Vi{net.fv.shape(0), 1}));

    JN_INFOA(net.o);

    JN_INFOS("cnn_ff done");
}

Arrayd cnn_bp(CNN &net, Arrayd &y) {
    Arrayd times({4}, 0);
    V<Timer> timers(4);

    int n = net.layers.size();
    auto &layer = net.layers.back();

    auto rot180 = [](const Arrayd &x) {
        int sz = x.size();
        Arrayd y(x.shape());

        const double *x_it = x.data() + sz - 1;
        double *y_it = y.data();

        for (int i = 0; i < sz; i++) {
            *y_it = *x_it;
            y_it++;
            x_it--;
        }
        return STD_ move(y);
        //return flipdim(flipdim(X, 1), 2);
    };

    JN_INFOS("cnn_bp start...");

    JN_INFOA(y);

    // error
    JN_INFOA(net.o);
    net.e = net.o - y;
    JN_INFOA(net.e);

    // loss function
    net.L = 1.0 / 2.0 * Arrayd::square(net.e).sum() / net.e.shape(1);
    JN_INFO(net.L);

    //  backprop deltas
    JN_INFOS("backprop deltas...");
    timers[0].toc();
    net.od = net.e * (net.o * (1 - net.o));   //  output delta
    JN_INFOA(net.od);
    JN_INFOS(string::format("%d-%d %d-%d", net.od.shape(0), net.od.shape(1), net.ffW.shape(0), net.ffW.shape(1)));
    net.fvd = multiply<double>(net.od, net.ffW, CblasNoTrans, CblasTrans);              //  feature vector delta
    JN_INFOA(net.fvd);
    if (layer.type == 'c') {         //  only conv layers has sigm function
        double *fvd_it = net.fvd.data();
        double *fv_it = net.fv.data();
        for (int i = 0; i != net.fvd.size(); i++) {
            *fvd_it *= (*fv_it) * (1 - (*fv_it));
            fvd_it++;
            fv_it++;
        }
        //net.fvd = net.fvd * (net.fv * (1 - net.fv));
        JN_INFOA(net.fvd);
    }
    times[0] += timers[0].toc();

    //  reshape feature vector deltas into output map style
    JN_INFOS("reshape feature vector deltas into output map style...");
    timers[1].toc();
    auto && sa = layer.a[0].shape();
    JN_INFOV(sa);
    int fvnum = sa[1] * sa[2];

    JN_INFO(fvnum);
    JN_INFO(layer.a.size());
    JN_INFOA(net.fvd);

    int asz = layer.a.size();

    // Set dj_it
    layer.d.resize(asz);
    V<double *> dj_it;
    for (int j = 0; j < asz; j++) {
        layer.d[j] = Arrayd(sa);
        dj_it.push_back(layer.d[j].data());
    }

    // Set dj
    double *fvd_it = net.fvd.data();
    for (int i = 0; i < sa[0]; i++) {
        for (int j = 0; j < asz; j++) {
            for (int k = 0; k < fvnum; k++) {
                *(dj_it[j]) = *fvd_it;
                fvd_it++;
                dj_it[j]++;
            }
        }
    }
    times[1] += timers[1].toc();

    JN_INFOS("l...");
    timers[2].toc();
    for (int l = n - 2; l >= 0; l--) {
        auto &p = net.layers[l];
        auto &q = net.layers[l+1];
        int scale = q.scale;
        int scale2 = scale*scale;
        if (p.type == 'c') {
            p.d.resize(p.a.size());
#ifdef USEMPI
            Shape shape = p.a[0].shape();
            auto mpialloc = mpi_alloc(p.a.size());
#pragma omp parallel for
            for (int j = mpialloc.beg; j < mpialloc.end; j++) {
#else
#pragma omp parallel for
            for (int j = 0; j < p.a.size(); j++) {
#endif
                p.d[j] = Arrayd(p.a[0].shape());
                auto && qd = expand(q.d[j], {1, q.scale, q.scale});
                double *d_it = p.d[j].data();
                double *a_it = p.a[j].data();
                double *qd_it = qd.data();
                for (int i = 0; i < qd.size(); i++) {
                    *d_it = (*a_it) * (1 - (*a_it)) * (*qd_it) / double(scale2);
                    d_it++;
                    a_it++;
                    qd_it++;
                }
                JN_INFO(j);
                JN_INFOA(p.d[j]);
            }
#ifdef USEMPI
            for (int j = 0; j < mpialloc.size; j++) {
                if (p.d[j].size() == 0) p.d[j] = Arrayd(shape);
                MPI_Bcast(p.d[j].data(), p.d[j].size(), MPI_DOUBLE, j/mpialloc.bin, MPI_COMM_WORLD);
            }
#endif
        }
        else if (p.type == 's') {
            p.d.resize(p.a.size());
#ifdef USEMPI
            Shape shape = p.a[0].shape();
            auto mpialloc = mpi_alloc(p.a.size());
#pragma omp parallel for
            for (int i = mpialloc.beg; i < mpialloc.end; i++) {
#else
#pragma omp parallel for
            for (int i = 0; i < p.a.size(); i++) {
#endif
                Arrayd z(p.a[0].shape(), 0);
                for (int j = 0; j < q.a.size(); j++) {
                    z.selfPlus(signal::convolve(q.d[j], rot180(q.k[i][j]), signal::CONVOLVE_FULL, signal::CONVOLVE_PLAIN));
                }
                p.d[i] = STD_ move(z);
                JN_INFO(i);
                JN_INFOA(p.d[i]);
            }
#ifdef USEMPI
            for (int j = 0; j < mpialloc.size; j++) {
                if (p.d[j].size() == 0) p.d[j] = Arrayd(shape);
                MPI_Bcast(p.d[j].data(), p.d[j].size(), MPI_DOUBLE, j/mpialloc.bin, MPI_COMM_WORLD);
            }
#endif
        }
    }
    times[2] += timers[2].toc();

    //  calc gradients
    JN_INFOS("calc gradients");
    timers[3].toc();
    for (int l = 1; l < n; l++) {
        auto &p = net.layers[l];
        auto &q = net.layers[l-1];
        if (p.type == 'c') {

            p.db.resize(p.a.size());
            p.dk.resize(q.a.size());
            JN_INFOS(string::format("p: %d, q: %d", p.a.size(), q.a.size()));

            V<Arrayd> ai;
            for (int i = 0; i < q.a.size(); i++) {
                p.dk[i].resize(p.a.size());
                ai.push_back(flipall(q.a[i]));
            }

#ifdef USEMPI
            Shape shape = signal::convolve(ai[0], p.d[0], signal::CONVOLVE_VALID, signal::CONVOLVE_PLAIN).shape();
            auto mpialloc = mpi_alloc(p.a.size());
#pragma omp parallel for
            for (int j = mpialloc.beg; j < mpialloc.end; j++) {
#else
#pragma omp parallel for
            for (int j = 0; j < p.a.size(); j++) {
#endif
                for (int i = 0; i < q.a.size(); i++) {
                    p.dk[i][j] = signal::convolve(ai[i], p.d[j], signal::CONVOLVE_VALID, signal::CONVOLVE_PLAIN);
                    p.dk[i][j].selfDivide(p.d[j].shape(0));
                    JN_INFOA(p.dk[i][j]);
                }
                p.db[j] = p.d[j].sum() / double(p.d[j].shape(0));
                JN_INFO(j);
                JN_INFO(p.db[j]);
            }
#ifdef USEMPI
            for (int j = 0; j < mpialloc.size; j++) {
                for (int i = 0; i < q.a.size(); i++) {
                    if (p.dk[i][j].size() == 0) p.dk[i][j] = Arrayd(shape);
                    MPI_Bcast(p.dk[i][j].data(), p.dk[i][j].size(), MPI_DOUBLE, j/mpialloc.bin, MPI_COMM_WORLD);
                }
                MPI_Bcast(&(p.db[j]), 1, MPI_DOUBLE, j/mpialloc.bin, MPI_COMM_WORLD);
            }
#endif
        }
    }
    times[3] += timers[3].toc();

    JN_INFOA(net.od);
    JN_INFOA(net.fv);

    net.dffW = multiply<double>(net.fv, net.od, CblasTrans, CblasNoTrans);
    net.dffW.selfDivide(double(net.od.shape(0)));
    net.dffb = meanv(net.od, 0);

    JN_INFOA(net.dffW);
    JN_INFOA(net.dffb);

    return times;

}

STD_ pair<double, Vi> cnn_test(CNN &net, Arrayd &x, Arrayd &y) {
    //  feedforward
    cnn_ff(net, x);

    JN_INFOA(net.o);
    JN_INFOA(y);

    Vi bad;
    for (int i = 0; i < y.shape(0); i++) {
        double max1 = net.o(i, 0);
        double max2 = y(i, 0);
        int max1_ind = 0, max2_ind = 0;
        for (int j = 0; j < y.shape(1); j++) {
            auto d = net.o(i, j);
            if (max1 < d) {
                max1 = d;
                max1_ind = j;
            }
            d = y(i, j);
            if (max2 < d) {
                max2 = d;
                max2_ind = j;
            }
        }
        if (max1_ind != max2_ind) bad.push_back(i);
    }

    double er = bad.size() / double(y.shape(0));

    JN_INFO(er);
    JN_INFOV(bad);

    return {er, bad};
}

void cnn_train(CNN &net, Arrayd &x, Arrayd &y, const CNNOpts &opts) {
    int i, l;

    int x_slis = x.shape(0);
    int x_rows = x.shape(1);
    int x_cols = x.shape(2);
    int x_area = x_rows*x_cols;

    int y_rows = y.shape(0);
    int y_cols = y.shape(1);

    int m = x.shape(0);
    if (m % opts.batchsize != 0) {
        ERROR_REPORT("Illegal numbatches!");
    }
    int numbatches = m / opts.batchsize;
    net.rL = {};
    Timer ff_timer, bp_timer, grad_timer, epoch_timer;
    Arrayd bp_times({4}, 0);
    double ff_time = 0, bp_time = 0, grad_time = 0;
    for (int i = 0; i < opts.numepochs; i++) {
        MPI_LOG << "[*] " << string::format("Start epoch %d/%d.", i, opts.numepochs) << std::endl;
        tic();

#ifdef RAND
        auto && kk = randperm(m);
#else
        auto && kk = randperm(m, m*100);
#endif
        JN_INFOA(kk);

        MPI_LOG << "[*] Batch" << std::endl;
        for (int l = 0; l < numbatches; l++) {
            MPI_LOG << " " << l << std::flush;
            JN_INFO(opts.batchsize);

            Arrayd batch_x({opts.batchsize, x_rows, x_cols});
            double *it_bx = batch_x.data();
            for (int i = 0; i < opts.batchsize; i++) {
                double *it_x = x.data() + kk[l*opts.batchsize+i]*x_area;
                for (int j = 0; j < x_area; j++) {
                    *it_bx = *it_x;
                    it_bx++;
                    it_x++;
                }
            }

            Arrayd batch_y({opts.batchsize, y_cols});
            double *it_by = batch_y.data();
            for (int i = 0; i < opts.batchsize; i++) {
                double *it_y = y.data() + kk[l*opts.batchsize+i]*y_cols;
                for (int j = 0; j < y_cols; j++) {
                    *it_by = *it_y;
                    it_by++;
                    it_y++;
                }
            }

            JN_INFOA(batch_x);
            JN_INFOA(batch_y);

            JN_INFOS("cnn_ff");
            ff_timer.toc();
            cnn_ff(net, batch_x);
            ff_time += ff_timer.toc();

            JN_INFOS("cnn_bp");
            bp_timer.toc();
            bp_times.selfPlus(cnn_bp(net, batch_y));
            bp_time += bp_timer.toc();

            JN_INFOS("cnn_applygrads");
            grad_timer.toc();
            cnn_applygrads(net, opts);
            grad_time += grad_timer.toc();

            if (net.rL.empty()) {
                net.rL.push_back(net.L);
            }
            net.rL.push_back(0.99 * net.rL.back() + 0.01 * net.L);
        }
        MPI_LOG << std::endl;
//        toc();
        JN_RUN(bp_times.print());
        MPI_LOG << "[*] " << string::format("epoch %d/%d: %f seconds.", i, opts.numepochs, epoch_timer.toc()) << std::endl;
        MPI_LOG << string::format("[*] ff time: %f, bp time: %f, grad time: %f seconds.", ff_time, bp_time, grad_time) << std::endl;
    }
}

void cnn_applygrads(CNN &net, const CNNOpts &opts) {
    int l, j, ii, sz;
    int incx = 1;
    int incy = 1;
    double a = double(-opts.alpha);
    for (int l = 1; l < net.layers.size(); l++) {
        auto &p = net.layers[l];
        auto &q = net.layers[l-1];
        if (p.type == 'c') {
#ifdef USEMPI
            auto mpialloc = mpi_alloc(p.a.size());
#pragma omp parallel for
            for (int j = mpialloc.beg; j < mpialloc.end; j++) {
#else
#pragma omp parallel for
            for (int j = 0; j < p.a.size(); j++) {
#endif
                for (int ii = 0; ii < q.a.size(); ii++) {
                    sz = p.k[ii][j].size();
                    JN_INFOA(p.dk[ii][j]);
                    daxpy(&sz, &a, p.dk[ii][j].data(), &incx, p.k[ii][j].data(), &incy);
                    JN_INFOA(p.k[ii][j]);
                    //p.k[ii][j] = p.k[ii][j] - opts.alpha * p.dk[ii][j];
                }
                p.b[j] = p.b[j] - opts.alpha * p.db[j];
            }
#ifdef USEMPI
            for (int j = 0; j < mpialloc.size; j++) {
                for (int i = 0; i < q.a.size(); i++) {
                    MPI_Bcast(p.dk[i][j].data(), p.dk[i][j].size(), MPI_DOUBLE, j/mpialloc.bin, MPI_COMM_WORLD);
                }
                MPI_Bcast(&(p.b[j]), 1, MPI_DOUBLE, j/mpialloc.bin, MPI_COMM_WORLD);
            }
#endif
        }
    }

    sz = net.ffW.size();
    for (int i = 0; i < sz; i++) net.ffW[i] -= opts.alpha * net.dffW[i];
    //net.ffW = net.ffW - opts.alpha * net.dffW;
    sz = net.ffb.size();
    for (int i = 0; i < sz; i++) net.ffb[i] -= opts.alpha * net.dffb[i];
    //net.ffb = net.ffb - opts.alpha * net.dffb;
    JN_INFOA(net.ffW);
    JN_INFOA(net.ffb);
}

void cnn_save(const CNN &cnn, Str filename) {
    STD_ ofstream ofile(filename.c_str());
    serialize(ofile, cnn);
    ofile.close();
}

Arrayd cnn_test_m(CNN &cnn, Arrayd &x) {
    cnn_ff(cnn, x);
    return cnn.o;
}

void serialize(STD_ ostream &stream, const CNNLayer &l) {
    ::cppsci::serialize(stream, l.type);
    ::cppsci::serialize(stream, l.a);
    ::cppsci::serialize(stream, l.d);
    ::cppsci::serialize(stream, l.b);
    ::cppsci::serialize(stream, l.db);
    ::cppsci::serialize(stream, l.k);
    ::cppsci::serialize(stream, l.dk);
    ::cppsci::serialize(stream, l.scale);
    ::cppsci::serialize(stream, l.outputmaps);
    ::cppsci::serialize(stream, l.kernelsize);

}

void parse(STD_ istream &stream, CNNLayer &l) {
    ::cppsci::parse(stream, l.type);
    ::cppsci::parse(stream, l.a);
    ::cppsci::parse(stream, l.d);
    ::cppsci::parse(stream, l.b);
    ::cppsci::parse(stream, l.db);
    ::cppsci::parse(stream, l.k);
    ::cppsci::parse(stream, l.dk);
    ::cppsci::parse(stream, l.scale);
    ::cppsci::parse(stream, l.outputmaps);
    ::cppsci::parse(stream, l.kernelsize);
}

void serialize(STD_ ostream &stream, const CNN &cnn) {
    ::cppsci::serialize(stream, cnn.layers);

    ::cppsci::serialize(stream, cnn.ffb);
    ::cppsci::serialize(stream, cnn.dffb);

    ::cppsci::serialize(stream, cnn.ffW);
    ::cppsci::serialize(stream, cnn.dffW);

    ::cppsci::serialize(stream, cnn.fv);
    ::cppsci::serialize(stream, cnn.fvd);

    ::cppsci::serialize(stream, cnn.e);
    ::cppsci::serialize(stream, cnn.o);
    ::cppsci::serialize(stream, cnn.od);

    ::cppsci::serialize(stream, cnn.rL);
    ::cppsci::serialize(stream, cnn.L);
}

void parse(STD_ istream &stream, CNN &cnn) {
    ::cppsci::parse(stream, cnn.layers);

    ::cppsci::parse(stream, cnn.ffb);
    ::cppsci::parse(stream, cnn.dffb);

    ::cppsci::parse(stream, cnn.ffW);
    ::cppsci::parse(stream, cnn.dffW);

    ::cppsci::parse(stream, cnn.fv);
    ::cppsci::parse(stream, cnn.fvd);

    ::cppsci::parse(stream, cnn.e);
    ::cppsci::parse(stream, cnn.o);
    ::cppsci::parse(stream, cnn.od);

    ::cppsci::parse(stream, cnn.rL);
    ::cppsci::parse(stream, cnn.L);
}

}} // namespace cppsci::CNN

