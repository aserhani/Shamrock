// -------------------------------------------------------//
//
// SHAMROCK code for hydrodynamics
// Copyright (c) 2021-2024 Timothée David--Cléris <tim.shamrock@proton.me>
// SPDX-License-Identifier: CeCILL Free Software License Agreement v2.1
// Shamrock is licensed under the CeCILL 2.1 License, see LICENSE for more information
//
// -------------------------------------------------------//

#pragma once

/**
 * @file PrimToConsGas.hpp
 * @author 
 * @brief
 * @date 2025-07-02
 */

#include "shambackends/vec.hpp"
#include "shamrock/solvergraph/IFieldSpan.hpp"
#include "shamrock/solvergraph/INode.hpp"
#include "shamrock/solvergraph/Indexes.hpp"

namespace shammodels::basegodunov::modules {
    template<class Tvec>
    class NodePrimToConsGas : public shamrock::solvergraph::INode {
        using Tscal = shambase::VecComponent<Tvec>;
        u32 block_size;
        Tscal gamma;

        public:
        NodePrimToConsGas(u32 block_size, Tscal gamma) : block_size(block_size), gamma(gamma) {}

        struct Edges {
            const shamrock::solvergraph::Indexes<u32> &sizes;
            const shamrock::solvergraph::IFieldSpan<Tvec> &spans_vel;
            const shamrock::solvergraph::IFieldSpan<Tscal> &spans_P;
            shamrock::solvergraph::IFieldSpan<Tscal> &spans_rho;
            shamrock::solvergraph::IFieldSpan<Tvec> &spans_rhov;
            shamrock::solvergraph::IFieldSpan<Tscal> &spans_rhoe;
        };

        inline void set_edges(
            std::shared_ptr<shamrock::solvergraph::Indexes<u32>> sizes,
            std::shared_ptr<shamrock::solvergraph::IFieldSpan<Tvec>> spans_vel,
            std::shared_ptr<shamrock::solvergraph::IFieldSpan<Tscal>> spans_P,
            std::shared_ptr<shamrock::solvergraph::IFieldSpan<Tscal>> spans_rho,
            std::shared_ptr<shamrock::solvergraph::IFieldSpan<Tvec>> spans_rhov,
            std::shared_ptr<shamrock::solvergraph::IFieldSpan<Tscal>> spans_rhoe) {
            __internal_set_ro_edges({sizes, spans_vel, spans_P});
            __internal_set_rw_edges({spans_rho, spans_rhov, spans_rhoe});
        }

        inline Edges get_edges() {
            return Edges{
                get_ro_edge<shamrock::solvergraph::Indexes<u32>>(0),
                get_ro_edge<shamrock::solvergraph::IFieldSpan<Tvec>>(1),
                get_ro_edge<shamrock::solvergraph::IFieldSpan<Tscal>>(2),
                get_rw_edge<shamrock::solvergraph::IFieldSpan<Tscal>>(0),
                get_rw_edge<shamrock::solvergraph::IFieldSpan<Tvec>>(1),
                get_rw_edge<shamrock::solvergraph::IFieldSpan<Tscal>>(2),
            };
        }

        void _impl_evaluate_internal();

        inline virtual std::string _impl_get_label() { return "PrimToConsGas"; };

        virtual std::string _impl_get_tex();
    };
} // namespace shammodels::basegodunov::modules
