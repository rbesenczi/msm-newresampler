/*
Copyright (c) 2022 King's College London, MeTrICS Lab, Renato Besenczi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
/*
 * Adaptive Barycentric code used with the permission of Tim Coulson under the below licence:
 * Copyright (C) 2014  Washington University School of Medicine
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the
 *  "Software"), to deal in the Software without restriction, including
 *  without limitation the rights to use, copy, modify, merge, publish,
 *  distribute, sublicense, and/or sell copies of the Software, and to
 *  permit persons to whom the Software is furnished to do so, subject to
 *  the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included
 *  in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 *  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 *  CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 *  TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 *  SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include "resampler.h"

using namespace std;

namespace newresampler {

MISCMATHS::FullBFMatrix Resampler::barycentric_data_interpolation(const Mesh& metric_in, const Mesh& sphLow){

    MISCMATHS::FullBFMatrix data(1, metric_in.nvertices());

    for (int i = 0; i < metric_in.nvertices(); i++)
        data.Set(1, i + 1, metric_in.get_pvalue(i));

    vector<map<int,double>> weights = get_adaptive_barycentric_weights(metric_in, sphLow);

    MISCMATHS::FullBFMatrix interpolated_data (1, sphLow.nvertices());

    for (int k = 0; k < sphLow.nvertices(); k++)
    {
        double val = 0.0;

        for (const auto& it: weights[k])
                val += data.Peek(1, it.first + 1) * it.second;

        interpolated_data.Set(1, k + 1, val);
    }

    return interpolated_data;
}

vector<std::map<int,double>> Resampler::get_adaptive_barycentric_weights(const Mesh& in_mesh, const Mesh& sphLow){
/*  ADAPTIVE BARYCENTRIC CODE SUPPLIED BY TIM COULSON AT WASHU Copyright (C) 2014  Washington University School of Medicine */
/* ORIGINAL CODE: https://github.com/Washington-University/workbench/blob/master/src/Files/SurfaceResamplingHelper.cxx      */

    Octree octreeSearch_in(in_mesh);
    vector<map<int,double>> forward = get_barycentric_weights(sphLow, in_mesh, octreeSearch_in);

    Octree octreeSearch_ref(sphLow);
    vector<map<int,double>> reverse = get_barycentric_weights(in_mesh, sphLow, octreeSearch_ref);

    vector<map<int,double>> reverse_reorder, adapt_weights;

    int numOldNodes = in_mesh.nvertices();
    int numNewNodes = sphLow.nvertices();
    vector<double> newAreas(numNewNodes, 0.0);
    vector<double> oldAreas(numOldNodes, 0.0);
    vector<double> correction(numOldNodes, 0.0);

    reverse_reorder.resize(forward.size());
    adapt_weights.resize(forward.size());

    for (int oldNode = 0; oldNode < numOldNodes; ++oldNode)//this loop can't be parallelized
    {
        oldAreas[oldNode] = compute_vertex_area(oldNode, in_mesh);
        for (auto iter = reverse[oldNode].begin();
             iter != reverse[oldNode].end(); ++iter)
            reverse_reorder[iter->first][oldNode] = iter->second;
    }

    for (int newNode = 0; newNode < numNewNodes; ++newNode)
    {
            newAreas[newNode] = compute_vertex_area(newNode, sphLow);
            if (reverse_reorder[newNode].size() <= forward[newNode].size())
                // choose whichever has the most samples as this will be the higher res mesh (this differs from source code)
                adapt_weights[newNode] = forward[newNode];
            else
                adapt_weights[newNode] = reverse_reorder[newNode];

            for (auto iter = adapt_weights[newNode].begin();
                 iter != adapt_weights[newNode].end(); ++iter)
            {
                //begin area correction by multiplying by target node calc_area
                iter->second *= newAreas[newNode];//begin the process of calc_area correction by multiplying by gathering node areas
                correction[iter->first] += iter->second;//now, sum the scattering weights to prepare for first normalization
            }
    }

    for (int newNode = 0; newNode < numNewNodes; ++newNode)
    {
        double weight_sum = 0.0;
            for (auto &iter: adapt_weights[newNode])//begin area correction by multiplying by target node calc_area
            {
                iter.second *= oldAreas[iter.first] / correction[iter.first];
                //divide the weights by their scatter sum, then multiply by current areas
                weight_sum += iter.second;//and compute the sum
            }
            if (weight_sum != 0.0)
                //this shouldn't happen unless no nodes remain due to roi, or node areas can be zero
                for (auto &iter: adapt_weights[newNode])
                    iter.second /= weight_sum;//and normalize to a sum of 1
    }

    return adapt_weights;
}

vector<std::map<int,double>> Resampler::get_barycentric_weights(const Mesh& low, const Mesh& orig, const Octree& oct){

    std::vector<std::map<int,double>> weights;

    for (int k = 0; k < low.nvertices(); k++)
    {
        Point v0, v1, v2;
        int n0, n1, n2;

        Point ci = low.get_coord(k);

        Triangle closest = oct.get_closest_triangle(ci);

        n0 = closest.get_vertex_no(0);
        n1 = closest.get_vertex_no(1);
        n2 = closest.get_vertex_no(2);
        v0 = closest.get_vertex_coord(0);
        v1 = closest.get_vertex_coord(1);
        v2 = closest.get_vertex_coord(2);

        weights.emplace_back(calc_barycentric_weights(v0, v1, v2, ci, n0, n1, n2));
    }

    return weights;
}

Mesh surface_resample(const Mesh& anatOrig, const Mesh& sphOrig, const Mesh& sphLow) {
//---ANATOMICAL MESH RESAMPLING---//
    Resampler resampler(Method::BARY);
    Mesh anatLow = sphLow;
    Octree octreeSearch(sphOrig);

    std::vector<map<int,double>> weights = resampler.get_barycentric_weights(sphLow, sphOrig, octreeSearch);

    for(int i = 0; i < sphLow.nvertices(); i++)
    {
        Point newPt;
        for (const auto& it: weights[i])
            newPt += anatOrig.get_coord(it.first) * it.second;
        anatLow.set_coord(i, newPt);
    }

    return anatLow;
}

Mesh metric_resample(const Mesh& metric_in, const Mesh& sphLow) {
//---METRIC FILE RESAMPLING---//

    Resampler resampler(Method::ADAP_BARY);

    MISCMATHS::FullBFMatrix interpolated_data = resampler.barycentric_data_interpolation(metric_in, sphLow);

    Mesh output = sphLow;

    for (int i = 0; i < sphLow.nvertices(); i++)
        output.set_pvalue(i, interpolated_data.Peek(1, i + 1));

    return output;
}

} //namespace newresampler
