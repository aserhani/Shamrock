// -------------------------------------------------------//
//
// SHAMROCK code for hydrodynamics
// Copyright (c) 2021-2024 Timothée David--Cléris <tim.shamrock@proton.me>
// SPDX-License-Identifier: CeCILL Free Software License Agreement v2.1
// Shamrock is licensed under the CeCILL 2.1 License, see LICENSE for more information
//
// -------------------------------------------------------//

/**
 * @file PrimToConsGas.cpp
 * @author 
 * @brief
 * @date 2025-07-02
 */

#include "shammodels/ramses/modules/PrimToConsGas.hpp"
#include "shambackends/kernel_call_distrib.hpp"
#include "shammath/riemann.hpp"
#include "shamrock/patch/PatchDataField.hpp"
#include "shamsys/NodeInstance.hpp"

namespace {

    template<class Tvec>
    struct KernelPrimToConsGas {
        using Tscal = shambase::VecComponent<Tvec>;

        inline static void kernel(
            const shambase::DistributedData<shamrock::PatchDataFieldSpanPointer<Tvec>> &spans_vel,
            const shambase::DistributedData<shamrock::PatchDataFieldSpanPointer<Tscal>> &spans_P,

            shambase::DistributedData<shamrock::PatchDataFieldSpanPointer<Tscal>> &spans_rho,
            shambase::DistributedData<shamrock::PatchDataFieldSpanPointer<Tvec>> &spans_rhov,
            shambase::DistributedData<shamrock::PatchDataFieldSpanPointer<Tscal>> &spans_rhoe,
            const shambase::DistributedData<u32> &sizes,
            u32 block_size,
            Tscal gamma) {

            shambase::DistributedData<u32> cell_counts
                = sizes.map<u32>([&](u64 id, u32 block_count) {
                      u32 cell_count = block_count * block_size;
                      return cell_count;
                  });

            sham::distributed_data_kernel_call(
                shamsys::instance::get_compute_scheduler_ptr(),
                sham::DDMultiRef{spans_vel, spans_P},
                sham::DDMultiRef{spans_rho, spans_rhov, spans_rhoe},
                cell_counts,
                [gamma](
                    u32 i,
                    const Tvec *__restrict vel,
                    const Tscal *__restrict P,
                    Tscal *__restrict rho,
                    Tvec *__restrict rhov,
                    Tscal *__restrict rhoe) {
                    auto prim_state = shammath::PrimState<Tvec>{rho[i], P[i], vel[i]};
                    auto cons_state = shammath::prim_to_cons(prim_state, gamma);

                    rho[i]   = cons_state.rho;
                    rhov[i]  = cons_state.rhovel;
                    rhoe[i]  = cons_state.rhoe;
                });
        }
    };

}

namespace shammodels::basegodunov::modules {

    template<class Tvec>
    void NodePrimToConsGas<Tvec>::_impl_evaluate_internal() {
        auto edges = get_edges();

        edges.spans_vel.check_sizes(edges.sizes.indexes); 
        edges.spans_P.check_sizes(edges.sizes.indexes);

        edges.spans_rho.ensure_sizes(edges.sizes.indexes);
        edges.spans_rhov.ensure_sizes(edges.sizes.indexes);
        edges.spans_rhoe.ensure_sizes(edges.sizes.indexes);

        KernelPrimToConsGas<Tvec>::kernel(
            edges.spans_vel.get_spans(),
            edges.spans_P.get_spans(),
            edges.spans_rho.get_spans(),
            edges.spans_rhov.get_spans(),
            edges.spans_rhoe.get_spans(),
            edges.sizes.indexes,
            block_size,
            gamma);
    }

    template<class Tvec>
    std::string NodePrimToConsGas<Tvec>::_impl_get_tex() {

        auto block_count = get_ro_edge_base(0).get_tex_symbol();
        auto vel = get_ro_edge_base(1).get_tex_symbol();
        auto P = get_ro_edge_base(2).get_tex_symbol();
        auto rho = get_rw_edge_base(0).get_tex_symbol();
        auto rhov = get_rw_edge_base(1).get_tex_symbol();
        auto rhoe = get_rw_edge_base(2).get_tex_symbol();

        std::string tex = R"tex(
            Primitive to conservative gas variables conversion

            \begin{align}
            % Todo
            \end{align}
        )tex";

        shambase::replace_all(tex, "{vel}", vel);
        shambase::replace_all(tex, "{P}", P);
        shambase::replace_all(tex, "{rho}", rho);
        shambase::replace_all(tex, "{rhov}", rhov);
        shambase::replace_all(tex, "{rhoe}", rhoe);
        shambase::replace_all(tex, "{block_count}", block_count);
        shambase::replace_all(tex, "{gamma}", shambase::format("{}", gamma));
        shambase::replace_all(tex, "{block_size}", shambase::format("{}", block_size));

        return tex;
    }

} // namespace shammodels::basegodunov::modules

template class shammodels::basegodunov::modules::NodePrimToConsGas<f64_3>;
