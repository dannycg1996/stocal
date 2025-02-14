"""Unit testing for rules in dsd.py """
import unittest
from stocal.tests.test_transitions import TestReactionRule as TestTransitionRule, TestMassAction


class TestBindingRule(unittest.TestCase):
    from stocal.examples.dsd import BindingRule
    Rule = BindingRule

    def test_lakin_r_b_example(self):
        # Test that the basic RB example from the Lakin paper can be replicated with the Binding Rule.
        r_b_1 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L' N^* R'}", "<L N^ R>")))[0].products.keys())[0]
        self.assertEqual(r_b_1, "{L'}<L>[N^]<R>{R'}")

    def test_lakin_r_b_example_diff_order(self):
        # Test that the basic RB example from the Lakin paper can be replicated with the Binding Rule regardless of input order.
        r_b_2 = list(list(set(self.Rule.novel_reactions(self.Rule(), "<L N^ R>", "{L' N^* R'}")))[0].products.keys())[0]
        self.assertEqual(r_b_2, "{L'}<L>[N^]<R>{R'}")

    def test_systems_which_can_bind_in_multiple_spots(self):
        # Tests that when possible, the Binding Rule yields multiple different bindings from the same inputs.
        r_b_3 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{S' N^* L' R'}", "<L N^ M N^>")))[0].products.keys())
        exp_res_3 = {"{S'}<L N^ M>[N^]{L' R'}", "{S'}<L>[N^]<M N^>{L' R'}"}
        self.assertEqual(set(), set.difference(r_b_3, exp_res_3))

    def test_binding_between_strands_where_the_output_has_no_lower_strand_before_the_double_strand(self):
        # Test a variant of the Binding Rule, where the yielded result doesn't have a lower strand preceding the d_s.
        r_b_4 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{ N^* L' R'}", "<L N^ M>")))[0].products.keys())[0]
        self.assertEqual(r_b_4, "<L>[N^]<M>{L' R'}")

    def test_binding_between_strands_where_the_output_has_no_lower_strand_after_the_double_strand(self):
        # Test a variant of the Binding Rule, where the yielded result doesn't have a lower strand after the d_s.
        r_b_5 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L' N^*}", "<L N^ M>")))[0].products.keys())[0]
        self.assertEqual(r_b_5, "{L'}<L>[N^]<M>")

    def test_binding_between_strands_where_the_output_has_no_upper_strand_before_the_double_strand(self):
        # Test a variant of the Binding Rule, where the yielded result doesn't have an upper strand preceding the d_s.
        r_b_6 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{A N^* L' R'}", "<N^ M>")))[0].products.keys())[0]
        self.assertEqual(r_b_6, "{A}[N^]<M>{L' R'}")

    def test_binding_between_strands_where_the_output_has_no_upper_strand_after_the_double_strand(self):
        # Test a variant of the Binding Rule, where the yielded result doesn't have an upper strand after the d_s.
        r_b_7 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L' N^* R}", "<L N^>")))[0].products.keys())[0]
        self.assertEqual(r_b_7, "{L'}<L>[N^]{R}")

    def test_simplest_binding_case(self):
        # Test the simplest strand to strand binding case, where the yielded result has just a single double toehold.
        r_b_8 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{N^*}", "<N^>")))[0].products.keys())[0]
        self.assertEqual(r_b_8, "[N^]")

    def test_lakin_fig_4a_example(self):
        # Test an example from Figure 4 of the Lakin paper
        r_b_9 = list(list(set(self.Rule.novel_reactions(self.Rule(), "<t^ x y>", "{t^*}[x]:[y u^]")))[0].products.keys())[0]
        self.assertEqual(r_b_9, "[t^]<x y>:[x]:[y u^]")

    def test_lakin_r_p_example(self):
        # Test that the basic RP example from the Lakin paper yields the correct result.
        r_b_10 = list(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 N^ S R1>", "{L' N^*}<L>[S R2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_b_10, "{L'}<L1>[N^]<S R1>:<L>[S R2]<R>{R'}")

    def test_binding_gate_to_gate_yields_no_results(self):
        # Test that binding does not occur between two gates.
        r_b_11 = set(self.Rule.novel_reactions(self.Rule(), "{N^* S' N^*}[C^]", "{L'}<L>[N^]<R>[M^]<S'>[A^]{B}"))
        self.assertEqual(r_b_11, set())

    def test_lower_strand_binding_to_gate(self):
        # Test that binding can occur between a lower strand and a gate.
        r_b_12 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{A C^*}", "{F}<B C^ G>[H^]<I>{J}")))[0].products.keys())[0]
        self.assertEqual(r_b_12, "{A}<B>[C^]::{F}<G>[H^]<I>{J}")

    def test_lower_strand_binding_to_second_gate(self):
        # Test that binding can occur between a lower strand and a gate, when the gate being bound to is preceded by another gate.
        r_b_13 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{F}<B C^ D G>[H^]:{J K}<I L>[M^]<N>{O}", "{A C^* E}")))[0].products.keys())[0]
        self.assertEqual(r_b_13, "{A}<B>[C^]{E}::{F}<D G>[H^]:{J K}<I L>[M^]<N>{O}")

    def test_upper_strand_binding_to_gate(self):
        # Test that binding can occur between an upper strand and a gate.
        r_b_14 = list(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 N^ S R1>", "{L' N^*}<L>[S R2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_b_14, "{L'}<L1>[N^]<S R1>:<L>[S R2]<R>{R'}")


class TestUnbindingRule(TestTransitionRule):
    from stocal.examples.dsd import UnbindingRule
    Rule = UnbindingRule

    def test_lakin_r_u_example(self):
        # r_u_1 tests that the basic RU example from the Lakin paper yields the correct result.
        r_u_1 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[N^]<R>{R'}")))[0].products.keys())
        exp_res_1 = {"{L' N^* R'}", "<L N^ R>"}
        self.assertEqual(set(), set.difference(r_u_1, exp_res_1))

    def test_unbinding_on_a_gate_containing_more_domains(self):
        # Test that RU correctly unbinds a gate which has more domains on its strands.
        r_u_2 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{B}<A>[D^]<C^ F>{C^* G}")))[0].products.keys())
        exp_res_2 = {"<A D^ C^ F>", "{B D^* C^* G}"}
        self.assertEqual(set(), set.difference(r_u_2, exp_res_2))

    def test_the_unbinding_of_the_second_gate_in_a_system(self):
        # Test a system which consists of two gates, with one possible point of unbinding, on the 2nd gate.
        r_u_3 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L1>[N^]<S R1>:<L>[S R2]<R>{R'}")))[0].products.keys())
        exp_res_3 = {"<L1 N^ S R1>", "{L' N^*}<L>[S R2]<R>{R'}"}
        self.assertEqual(set(), set.difference(r_u_3, exp_res_3))

    def test_the_unbinding_of_a_system_with_several_possible_unbinding_locations(self):
        # Test a system which can unbind at 3 different points.
        r_u_4 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{A}<B>[C^]<D>{E}::{F}<G>[H^]<I>{J}::{K}<L>[M^]<N>{O}")))[0].products.keys())
        exp_res_4 = {"{F}<B C^ D G>[H^]{J}::{K}<I L>[M^]<N>{O}", "{A C^* E}", "{A}<B>[C^]{E}::{K}<D G H^ I L>[M^]<N>{O}", "{F H^* J}",
                     "{A}<B>[C^]{E}::{F}<D G>[H^]<I L M^ N>{J}", "{K M^* O}"}
        self.assertEqual(set(), set.difference(r_u_4, exp_res_4))


class TestCoveringRule(TestTransitionRule):
    from stocal.examples.dsd import CoveringRule
    Rule = CoveringRule

    def test_lakin_r_c_example_l_to_r(self):
        # Tests that the basic RC  example from the Lakin paper yields the correct result.
        r_c_1 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S]<N^ R>{N^* R'}")))[0].products.keys())[0]
        self.assertEqual(r_c_1, "{L'}<L>[S N^]<R>{R'}")

    def test_lakin_rc_example_r_to_l(self):
        # r_c_2 tests that the RC example works in reverse, in the right to left direction.
        r_c_2 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L' N^*}<L N^>[S]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_c_2, "{L'}<L>[N^ S]<R>{R'}")

    def test_covering_rule_variant_left_to_right(self):
        # Test a basic variant of the covering rule RC, applied left to right.
        r_c_3 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[S]<N^ R>{N^* R'}")))[0].products.keys())[0]
        self.assertEqual(r_c_3, "[S N^]<R>{R'}")

    def test_covering_rule_variant_right_to_left(self):
        # Test a basic variant of the covering rule RC, applied right to left.
        r_c_4 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{R' N^*}<R N^>[S]")))[0].products.keys())[0]
        self.assertEqual(r_c_4, "{R'}<R>[N^ S]")

    def test_covering_rule_across_gates_which_are_joined_via_upper_strand(self):
        # Test the application of the covering rule across gates, left to right, where the gates are joined by an upper strand.
        r_c_5 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{A}<B>[C]{E^*}::{F}<E^ D>[G]")))[0].products.keys())[0]
        self.assertEqual(r_c_5, "{A}<B>[C E^]::{F}<D>[G]")
        # N.B: No right_to_left version of this exists, due to the chosen normal form.

    def test_covering_rule_across_gates_which_are_joined_via_upper_strand_variant(self):
        # A variation of the last test, where the lower domain which is being bound to is followed by other domains.
        r_c_6 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{A}<B>[C]{E^* Z}::{F}<E^ D>[G]")))[0].products.keys())[0]
        self.assertEqual(r_c_6, "{A}<B>[C E^]{Z}::{F}<D>[G]")
        # N.B: No right_to_left version of this exists, due to the chosen normal form.

    def test_covering_rule_left_to_right_variant(self):
        # Tests a variation of the covering rule where the gate which is being 'covered' is followed immediately by another d_s.
        r_c_7 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S]<N^ R>{N^* R'}::[A B]")))[0].products.keys())[0]
        self.assertEqual(r_c_7, "{L'}<L>[S N^]<R>{R'}::[A B]")
        # N.B: No right_to_left version of this exists, due to the chosen normal form.

    def test_covering_rule_left_to_right_variant_2(self):
        # Tests a variation of the covering rule where the gate which is being 'covered' lies between other gates.
        r_c_8 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[C D]<A>:{L'}<L>[S]<N^ R>{N^* R'}::[A B]")))[0].products.keys())[0]
        self.assertEqual(r_c_8, "[C D]<A>:{L'}<L>[S N^]<R>{R'}::[A B]")
        # N.B: No right_to_left version of this exists, due to the chosen normal form.


class TestMigrationRule(TestTransitionRule):
    from stocal.examples.dsd import MigrationRule
    Rule = MigrationRule

    def test_lakin_r_m_example_upper_l_to_r(self):
        # r_m_1 tests that the basic RM example from the Lakin paper yields the correct result.
        r_m_1 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S R2>:<L1>[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_m_1, "{L'}<L>[S1 S]<R2>:<L1 S>[S2]<R>{R'}")

    def test_lakin_r_m_example_lower_l_to_r(self):
        # Test variants of r_m_1 but when the overhang is on the lower strand:
        r_m_2 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S R2}::{L1}[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_m_2, "{L'}<L>[S1 S]{R2}::{L1 S}[S2]<R>{R'}")

    def test_lakin_r_m_example_upper_r_to_l(self):
        # Tests that the basic RM example from the Lakin paper yields the correct result - when done in reverse (right to left).
        r_m_3 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1 S]<R2>:<L1 S>[S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_m_3, "{L'}<L>[S1]<S R2>:<L1>[S S2]<R>{R'}")

    def test_lakin_r_m_example_lower_r_to_l(self):
        # Tests that the lower strand version of the RM example from the Lakin paper can be performed left to right (reverse)
        r_m_4 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1 S]<R2>:<L1 S>[S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_m_4, "{L'}<L>[S1]<S R2>:<L1>[S S2]<R>{R'}")

    def test_lakin_r_m_example_upper_l_to_r_second_overhang_only_in_result(self):
        # Test variant of r_m_1 where R2 is missing (so the result only has one overhang):
        r_m_5 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S>:<L1>[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_m_5, "{L'}<L>[S1 S]:<L1 S>[S2]<R>{R'}")

    def test_lakin_r_m_example_upper_r_to_l_second_overhang_only_in_input(self):
        # Test variant of RM (applied right to left) where the input only has the 2nd overhang. Also reverse of r_m_5.
        r_m_6 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1 S]:<L1 S>[S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_m_6, "{L'}<L>[S1]<S>:<L1>[S S2]<R>{R'}")

    def test_lakin_r_m_example_lower_l_to_r_second_overhang_only_in_result(self):
        # Test variant of r_m_2 where R2 is missing (so the result only has one overhang):
        r_m_7 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S}::{L1}[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_m_7, "{L'}<L>[S1 S]::{L1 S}[S2]<R>{R'}")

    def test_lakin_r_m_example_lower_r_to_l_second_overhang_only_in_input(self):
        # Test lower strand variant of RM (applied right to left) where the input only has the 2nd overhang. Also reverse of r_m_7.
        r_m_8 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1 S]::{L1 S}[S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_m_8, "{L'}<L>[S1]{S}::{L1}[S S2]<R>{R'}")

    def test_lakin_r_m_example_upper_l_to_r_input_only_has_first_overhang(self):
        # Test variant of r_m_1 where the input only has the 1st overhang (i.e. L1 is missing)
        r_m_9 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S R2>:[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_m_9, "{L'}<L>[S1 S]<R2>:<S>[S2]<R>{R'}")

    def test_lakin_r_m_example_upper_r_to_l_result_only_has_first_overhang(self):
        # Test r_m_9 applied in reverse (right to left) where the result only has the 1st overhang.
        r_m_10 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1 S]<R2>:<S>[S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_m_10, "{L'}<L>[S1]<S R2>:[S S2]<R>{R'}")

    def test_lakin_r_m_example_lower_l_to_r_input_only_has_first_overhang(self):
        # Test lower strand variant of r_m_1 where the input only has the 1st overhang (i.e. L1 is missing).
        r_m_11 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S R2}::[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_m_11, "{L'}<L>[S1 S]{R2}::{S}[S2]<R>{R'}")

    def test_lakin_r_m_example_lower_r_to_l_result_only_has_first_overhang(self):
        # Test lower strand variant of RM (appied right to left) where the result only has the 1st overhang. Also reverse of r_m_11
        r_m_12 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1 S]{R2}::{S}[S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_m_12, "{L'}<L>[S1]{S R2}::[S S2]<R>{R'}")

    def test_lakin_r_m_example_upper_l_to_r_input_only_has_the_first_overhang_and_result_only_has_second_overhang(self):
        # Test variants of r_m_1 where R2 and L1 are missing:
        r_m_13 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S>:[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_m_13, "{L'}<L>[S1 S]:<S>[S2]<R>{R'}")

    def test_lakin_r_m_example_upper_r_to_l_input_only_has_the_second_overhang_and_result_only_has_first_overhang(self):
        # Test variant of Lakin's RM rule (applied right to left) where R2 and L1 are missing. Also reverse of r_m_13.
        r_m_14 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1 S]:<S>[S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_m_14, "{L'}<L>[S1]<S>:[S S2]<R>{R'}")

    def test_lakin_r_m_example_lower_l_to_r_input_only_has_the_first_overhang_and_result_only_has_second_overhang(self):
        # Test variants of r_m_2 where R2 and L1 are missing:
        r_m_15 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S}::[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_m_15, "{L'}<L>[S1 S]::{S}[S2]<R>{R'}")

    def test_lakin_r_m_example_lower_r_to_l_input_only_has_the_second_overhang_and_result_only_has_first_overhang(self):
        # Test lower strand variant of Lakin's RM rule (applied right to left) where R2 and L1 are missing. Also reverse of r_m_15.
        r_m_16 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1 S]::{S}[S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(r_m_16, "{L'}<L>[S1]{S}::[S S2]<R>{R'}")

    def test_that_migration_rule_is_not_applied_to_lakin_displacement_example_rd(self):
        # Test that RM is not applied on the RD example, as the two should be mutually exclusive.
        r_m_17 = set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S R>:<L2>[S]<R2>{R'}"))
        self.assertEqual(r_m_17, set())

    def test_that_migration_rule_is_not_applied_to_lower_strand_version_of_lakin_displacement_example_rd(self):
        # Test that RM is not applied to the lower strand version of the RD example, as the rules should be mutually exclusive.
        r_m_18 = set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S R}::{L2}[S]<R2>{R'}"))
        self.assertEqual(r_m_18, set())

    def test_that_migration_rule_is_not_applied_to_lakin_displacement_example_fig_4a(self):
        # Test that the RM rule is not applied to the RD example from Figure 4a).
        r_m_19 = set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x]:[y u^]"))
        self.assertEqual(r_m_19, set())

    def test_that_migration_rule_is_not_applied_to_lower_strand_version_of_lakin_displacement_example_fig_4a(self):
        # Test that the RM rule is not applied to the lower strand version of the RD example from Figure 4a).
        r_m_20 = set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::[x]::[y u^]"))
        self.assertEqual(r_m_20, set())

    def test_upper_l_to_r_lakin_fig_4a_migration_example_correct(self):
        # Test the migration rule is applied correctly to the example from Figure 4a) of Lakin's paper.
        r_m_21 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]<y>:[y u^]")))[0].products.keys())[0]
        self.assertEqual(r_m_21, "[t^ x y]:<y>[u^]")

    def test_upper_r_to_l_lakin_fig_4a_migration_example_correct(self):
        # Test the migration rule is applied correctly (in reverse) to the example from Figure 4a) of Lakin's paper (i.e. r_m_21).
        r_m_22 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x y]:<y>[u^]")))[0].products.keys())[0]
        self.assertEqual(r_m_22, "[t^ x]<y>:[y u^]")

    def test_lower_l_to_r_lakin_fig_4a_migration_example_correct(self):
        # Test the migration rule is applied correctly to the lower strand version of the example from Figure 4a) of Lakin's paper.
        r_m_23 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]{y}::[y u^]")))[0].products.keys())[0]
        self.assertEqual(r_m_23, "[t^ x y]::{y}[u^]")

    def test_lower_r_to_l_lakin_fig_4a_migration_example_correct(self):
        # Test that the rule works (right-to-left) on the lower strand version of the Fig. 4a example (i.e. r_m_23) in Lakin's paper.
        r_m_24 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x y]::{y}[u^]")))[0].products.keys())[0]
        self.assertEqual(r_m_24, "[t^ x]{y}::[y u^]")

    def test_migration_rule_upper_l_to_r_variant_1(self):
        # Test system where the 2nd gate involved in migration is connected to a 3rd gate via the upper strand.
        r_m_25 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x v]::[y u^]")))[0].products.keys())[0]
        self.assertEqual("[t^ x]<y>:<x>[v]::[y u^]", r_m_25)

    def test_migration_rule_upper_r_to_l_variant_1(self):
        # Test right-to-left rule application where 2nd gate is connected to a 3rd via an upper strand. Reverse of r_m_25.
        r_m_26 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]<y>:<x>[v]::[y u^]")))[0].products.keys())[0]
        self.assertEqual("[t^]<x y>:[x v]::[y u^]", r_m_26)

    def test_migration_rule_lower_l_to_r_variant_1(self):
        # Test system where the 2nd gate involved in migration is connected to a 3rd gate via the lower strand.
        r_m_27 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::[x v]:[y u^]")))[0].products.keys())[0]
        self.assertEqual("[t^ x]{y}::{x}[v]:[y u^]", r_m_27)

    def test_migration_rule_lower_r_to_l_variant_1(self):
        # Test right-to-left rule application of a system where the 2nd gate involved connects to a 3rd gate via the lower strand.
        # Reverse of r_m_27
        r_m_28 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]{y}::{x}[v]:[y u^]")))[0].products.keys())[0]
        self.assertEqual("[t^]{x y}::[x v]:[y u^]", r_m_28)


class TestDisplacementRule(TestTransitionRule):
    from stocal.examples.dsd import DisplacementRule
    Rule = DisplacementRule

    def test_lakin_r_d_example_upper_l_to_r(self):
        # Test the rule reduction example RD from Lakin's paper.
        r_d_1 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S R>:<L2>[S]<R2>{R'}")))[0].products.keys())
        exp_res_1 = {"<L2 S R2>", "{L'}<L>[S1 S]<R>{R'}"}
        self.assertEqual(set(), set.difference(r_d_1, exp_res_1))

    def test_lakin_r_d_example_upper_r_to_l(self):
        # Test an inverted version of example RD (r_d_1 above) from Lakin's paper, where the rule is applied right to left.
        r_d_2 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S]<L2>:<R S>[S1]<R2>{R'}")))[0].products.keys())
        exp_res_2 = {"<L S L2>", "{L'}<R>[S S1]<R2>{R'}"}
        self.assertEqual(set(), set.difference(r_d_2, exp_res_2))

    def test_lakin_r_d_example_lower_l_to_r(self):
        # Test the lower strand equivalent of the reduction example RD (r_d_1 above) from Lakin's paper.
        r_d_3 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S R}::{L2}[S]<R2>{R'}")))[0].products.keys())
        exp_res_3 = {"{L2 S R'}", "{L'}<L>[S1 S]<R2>{R}"}
        self.assertEqual(set(), set.difference(r_d_3, exp_res_3))

    def test_lakin_r_d_example_lower_r_to_l(self):
        # Test an inverted lower strand version of example RD (r_d_1 above) from Lakin's paper, applying the rule right-to-left.
        r_d_4 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S]{L2}::{R S}[S1]<R2>{R'}")))[0].products.keys())
        exp_res_4 = {"{L' S L2}", "{R}<L>[S S1]<R2>{R'}"}
        self.assertEqual(set(), set.difference(r_d_4, exp_res_4))

    def test_lakin_fig_4a_example_upper_l_to_r(self):
        # Tests that the application of the displacement rule from Figure 4a works as expected.
        r_d_5 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x]:[y u^]")))[0].products.keys())
        exp_res_5 = {"<x>", "[t^ x]<y>:[y u^]"}
        self.assertEqual(set(), set.difference(r_d_5, exp_res_5))

    def test_lakin_fig_4a_example_upper_r_to_l(self):
        # Tests that an altered version of the displacement eg. from Fig 4a can be displaced in the right-to-left direction.
        r_d_6 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[u^ y]:[x]:<y x>[t^]")))[0].products.keys())
        exp_res_6 = {"<x>", "[u^ y]:<y>[x t^]"}
        self.assertEqual(set(), set.difference(r_d_6, exp_res_6))

    def test_lakin_fig_4a_example_lower_l_to_r(self):
        # Tests that the application of the Displacement example from Figure 4a works as expected.
        r_d_7 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::[x]::[y u^]")))[0].products.keys())
        exp_res_7 = {"{x}", "[t^ x]{y}::[y u^]"}
        self.assertEqual(set(), set.difference(r_d_7, exp_res_7))

    def test_lakin_fig_4a_example_lower_r_to_l(self):
        # Tests an inverted (lower strand) version of the displacement example from Fig 4a (in the right-to-left direction).
        r_d_8 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[u^ y]::[x]::{y x}[t^]")))[0].products.keys())
        exp_res_8 = {"{x}", "[u^ y]::{y}[x t^]"}
        self.assertEqual(set(), set.difference(r_d_8, exp_res_8))

    def test_lakin_migration_example_fig_upper_4a_l_to_r_does_not_yield_results(self):
        # Test that the Displacement rule does not get applied to the Migration example from Figure 4a of the Lakin paper.
        r_d_9 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]<y>:[y u^]"))))
        self.assertEqual(set(), r_d_9)

    def test_lakin_migration_example_fig_upper_4a_r_to_l_does_not_yield_results(self):
        # Tests that this rule yields no results when applied to an inverted Migration example from Fig. 4a of the Lakin paper.
        r_d_10 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[u^ y]:<y>[x t]"))))
        self.assertEqual(set(), r_d_10)

    def test_lakin_migration_example_fig_lower_4a_l_to_r_does_not_yield_results(self):
        # Test that the lower strand version of the example from Fig. 4a cannot yield displacement products.
        r_d_11 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]{y}::[y u^]"))))
        self.assertEqual(r_d_11, set())

    def test_lakin_migration_example_fig_lower_4a_r_to__does_not_yield_results(self):
        # Tests that this rule yields no results when applied to an inverted, flipped Migration example from Fig. 4a of the Lakin paper.
        r_d_12 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[u^ y]::{y}[x t]"))))
        self.assertEqual(set(), r_d_12)

    def test_that_more_migration_examples_yield_no_displacement_results(self):
        # Test that other systems where migration can occur cannot be displaced:
        r_d_13 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x v]::[y u^]"))))
        self.assertEqual(set(), r_d_13)

    def test_that_more_migration_examples_yield_no_displacement_results_2(self):
        # Test that other systems where migration can occur cannot be displaced:
        r_d_14 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::[x v]:[y u^]"))))
        self.assertEqual(set(), r_d_14)

    def test_displacement_of_upper_strand_which_connects_to_the_next_gate_via_upper_strand_l_to_r(self):
        # This test checks that applying the displacement rule along an upper strand works, when the strand which is being
        # displaced is connected along its upper strand to the next gate (left to right).
        r_d_15 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x]::[y u^]")))[0].products.keys())
        exp_res_15 = {"[t^ x]<y>", "<x>[y u^]"}
        self.assertEqual(set(), set.difference(r_d_15, exp_res_15))

    def test_displacement_of_upper_strand_which_connects_to_the_previous_gate_via_upper_strand_r_to_l(self):
        # This test checks that applying the displacement rule along an upper strand works, when the strand which is being
        # displaced is connected along its upper strand to the previous gate (right to left). Variant of r_d_15.
        r_d_16 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[u^ y]::[x]:<y x>[t^]")))[0].products.keys())
        exp_res_16 = {"[u^ y]<x>", "<y>[x t^]"}
        self.assertEqual(set(), set.difference(r_d_16, exp_res_16))

    def test_displacement_of_lower_strand_which_connects_to_the_next_gate_via_lower_strand_l_to_r(self):
        # This test checks that applying the displacement rule along a lower strand works, when the strand which is being
        # displaced is connected to the next gate (left to right) along its lower strand.
        r_d_17 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::[x]:[y u^]")))[0].products.keys())
        exp_res_17 = {"[t^ x]{y}", "{x}[y u^]"}
        self.assertEqual(set(), set.difference(r_d_17, exp_res_17))

    def test_displacement_of_lower_strand_which_connects_to_the_previous_gate_via_lower_strand_r_to_l(self):
        # This test checks that applying the displacement rule along an lower strand works, when the toehold which is being
        # displaced is connected along its upper strand to the previous gate (right to left). Variant of r_d_16
        r_d_18 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[u^ y]:[x]::{y x}[t^]")))[0].products.keys())
        exp_res_18 = {"[u^ y]{x}", "{y}[x t^]"}
        self.assertEqual(set(), set.difference(r_d_18, exp_res_18))

    def test_displacement_of_upper_strand_which_is_connected_to_the_next_strand_via_upper_strand_l_to_r_variant_1(self):
        # This tests that displacing an upper strand works, when the strand which is being displaced is connected along to the
        # next gate (left to right) via the upper strand. Variant of r_d_15 but with an upper strand attached to the second d_s.
        r_d_19 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:<R>[x]::[y u^]")))[0].products.keys())
        exp_res_19 = {"[t^ x]<y>", "<R x>[y u^]"}
        self.assertEqual(set(), set.difference(r_d_19, exp_res_19))

    def test_displacement_of_upper_strand_which_is_connected_to_the_previous_strand_via_upper_strand_r_to_l_variant_1(self):
        # This tests that displacing an upper strand (right to left) works, when the strand which is being displaced is connected
        # along to the previous gate via the upper strand. Variant of r_d_16 but with an upper strand attached to the second d_s.
        r_d_20 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[u^ y]::[x]<R>:<y x>[t^]")))[0].products.keys())
        exp_res_20 = {"[u^ y]<x R>", "<y>[x t^]"}
        self.assertEqual(set(), set.difference(r_d_20, exp_res_20))

    def test_displacement_of_lower_strand_which_is_connected_to_the_next_strand_via_lower_strand_l_to_r_variant_1(self):
        # This tests that displacing a lower strand works, when the strand which is being displaced is connected along to the
        # next gate (left to right) via a lower strand. Variant of r_d_17 but with a lower strand attached to the second d_s.
        r_d_21 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::{R}[x]:[y u^]")))[0].products.keys())
        exp_res_21 = {"[t^ x]{y}", "{R x}[y u^]"}
        self.assertEqual(set(), set.difference(r_d_21, exp_res_21))

    def test_displacement_of_lower_strand_which_is_connected_to_the_previous_strand_via_lower_strand_r_to_l_variant_1(self):
        # This tests that displacing a lower strand (right-to-left) works, when the strand which is being displaced is connected
        # to the previous gate via a lower strand. Variant of r_d_18 but with a lower strand attached to the second d_s.
        r_d_22 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[u^ y]:[x]{R}::{y x}[t^]")))[0].products.keys())
        exp_res_22 = {"[u^ y]{x R}", "{y}[x t^]"}
        self.assertEqual(set(), set.difference(r_d_22, exp_res_22))

    def test_displacement_of_upper_strand_which_is_connected_to_the_next_strand_via_upper_strand_l_to_r_variant_2(self):
        # This tests that displacing an upper strand (left-to-right) works, when the strand which is being displaced is connected
        # to the next gate via the upper strand. Variant of r_d_19 but with a lower strand attached to the second d_s.
        r_d_23 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:<r>[x]{g}::[y u^]")))[0].products.keys())
        exp_res_23 = {"[t^ x]<y>{g}", "<r x>[y u^]"}
        self.assertEqual(set(), set.difference(r_d_23, exp_res_23))

    def test_displacement_of_upper_strand_which_is_connected_to_the_previous_strand_via_upper_strand_r_to_l_variant_2(self):
        # This tests that displacing an upper strand (right-to-right) works, when the strand which is being displaced is connected
        # to the previous gate via the upper strand. Variant of r_d_20 but with a lower strand attached to the first d_s.
        r_d_24 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[u^ y]::{g}[x]<r>:<y x>[t^]")))[0].products.keys())
        exp_res_24 = {"[u^ y]<x r>", "{g}<y>[x t^]"}
        self.assertEqual(set(), set.difference(r_d_24, exp_res_24))

    def test_displacement_of_lower_strand_which_is_connected_to_the_next_strand_via_lower_strand_l_to_r_variant_2(self):
        # This tests that displacing a lower strand (left-to-right) works, when the strand which is being displaced is connected
        # to the next gate via the upper strand. Variant of r_d_21 but with a upper strand attached to the second d_s.
        r_d_25 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::{r}[x]<g>:[y u^]")))[0].products.keys())
        exp_res_25 = {"[t^ x]<g>{y}", "{r x}[y u^]"}
        self.assertEqual(set(), set.difference(r_d_25, exp_res_25))

    def test_displacement_of_lower_strand_which_is_connected_to_the_previous_strand_via_lower_strand_r_to_l_variant_2(self):
        # This tests that displacing a lower strand (right-to-left) works, when the strand which is being displaced is connected
        # to the previous gate via the lower strand. Variant of r_d_22 but with an upper strand attached to the first d_s.
        r_d_26 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[u^ y]:<g>[x]{r}::{y x}[t^]")))[0].products.keys())
        exp_res_26 = {"[u^ y]{x r}", "{y}<g>[x t^]"}
        self.assertEqual(set(), set.difference(r_d_26, exp_res_26))


class TestStrandLeakageRule(unittest.TestCase):
    from stocal.examples.dsd import StrandLeakageRule
    Rule = StrandLeakageRule

    def test_lakin_l_s_example(self):
        # Test that the basic LS example from the Lakin paper can be replicated with the Leakage Rule.
        l_s_1 = set(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 S R1>", "{L'}<L>[S]<R>{R'}")))[0].products.keys())
        exp_res_1 = {"<L S R>", "{L'}<L1>[S]<R1>{R'}"}
        self.assertEqual(set(), set.difference(l_s_1, exp_res_1))

    def test_lakin_l_s_example_rotated(self):
        # Test the basic LS example from the Lakin paper, but rotate the invader strand to be a lower strand.
        l_s_2 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L1 S* R1}", "{L'}<L>[S]<R>{R'}")))[0].products.keys())
        exp_res_2 = {"{L' S* R'}", "{L1}<L>[S]<R>{R1}"}
        self.assertEqual(set(), set.difference(l_s_2, exp_res_2))

    def test_that_strand_leakage_does_not_apply_to_short_double_toeholds(self):
        # Test that the strand leakage rule yields nothing when a gate's double strand has form [N^].
        l_s_3 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L1 S* R1}", "{L'}<L>[S^]<R>{R'}"))))
        self.assertEqual(set(), l_s_3)

    def test_that_strand_leakage_fails_when_invader_strand_does_not_match_gate(self):
        # Test that when the invader sequence of domains does not match the sequence of domains within the d_s of the
        # other input, no leakages are yielded
        l_s_4 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L1 A* B^* C* R1}", "{L'}<L>[A B C]<R>{R'}"))))
        self.assertEqual(set(), l_s_4)

    def test_strand_leakage_with_an_upper_invader_which_causes_a_gate_to_leak_its_upper_strand(self):
        # Test the LS rule when the invader strand is an upper strand which contains a mixture of toeholds and long domains.
        l_s_5 = set(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 S T^ R1>", "{L'}<L>[S T^]<R>{R'}")))[0].products.keys())
        exp_res_5 = {"<L S T^ R>", "{L'}<L1>[S T^]<R1>{R'}"}
        self.assertEqual(set(), set.difference(l_s_5, exp_res_5))

    def test_strand_leakage_with_an_upper_invader_which_causes_a_gate_to_leak_its_lower_strand(self):
        # Test the LS rule when the invader strand is an upper strand which can only initiate a leakage after rotating into
        # a lower strand.
        l_s_6 = set(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 T^* S* R1>", "{L'}<L>[S T^]<R>{R'}")))[0].products.keys())
        exp_res_6 = {"{L' S* T^* R'}", "{R1}<L>[S T^]<R>{L1}"}
        self.assertEqual(set(), set.difference(l_s_6, exp_res_6))

    def test_strand_leakage_with_a_lower_invader_which_causes_a_gate_to_leak_its_lower_strand(self):
        # Test the LS rule when the invader strand is a lower strand which contains a mixture of toeholds and long domains.
        l_s_7 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L1 S* T^* R1}", "{L'}<L>[S T^]<R>{R'}")))[0].products.keys())
        exp_res_7 = {"{L' S* T^* R'}", "{L1}<L>[S T^]<R>{R1}"}
        self.assertEqual(set(), set.difference(l_s_7, exp_res_7))

    def test_strand_leakage_with_a_lower_invader_which_causes_a_gate_to_leak_its_upper_strand(self):
        # Test the LS rule when the invader strand is a lower strand which can only initiate a leakage after rotating into
        # an upper strand.
        l_s_8 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L1 T^ S R1}", "{L'}<L>[S T^]<R>{R'}")))[0].products.keys())
        exp_res_8 = {"<L S T^ R>", "{L'}<R1>[S T^]<L1>{R'}"}
        self.assertEqual(set(), set.difference(l_s_8, exp_res_8))

    def test_strand_leakage_with_constructs_which_contain_more_complex_sequences_of_domains_1(self):
        # Test the LS rule with an upper invader strand which can only cause a leak with one rotation i.e. if the invader rotates
        # into a lower strand, a leakage will not occur (on the lower strand). Variant of l_s_5 with long sequences of domains.
        l_s_9 = set(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 LA S T^ RA R1>",
                    "{L' L2}<L LB>[S T^]<RB R>{R2 R'}")))[0].products.keys())
        exp_res_9 = {"<L LB S T^ RB R>", "{L' L2}<L1 LA>[S T^]<RA R1>{R2 R'}"}
        self.assertEqual(set(), set.difference(l_s_9, exp_res_9))

    def test_strand_leakage_with_constructs_which_contain_more_complex_sequences_of_domains_1(self):
        # Test the LS rule when the invader strand is an upper strand which can only cause a leak if it rotates into a lower strand.
        # Variant of l_s_6 with longer sequences of domains.
        l_s_10 = set(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 LA T^* S* RA R1>",
                                                        "{L' L2}<L LB>[S T^]<RB R>{R2 R'}")))[0].products.keys())
        exp_res_10 = {"{L' L2 S* T^* R2 R'}", "{R1 RA}<L LB>[S T^]<RB R>{LA L1}"}
        self.assertEqual(set(), set.difference(l_s_10, exp_res_10))

    def test_strand_leakage_with_constructs_which_contain_more_complex_sequences_of_domains_2(self):
        # Test the LS rule when the invader is a lower strand which can only cause a leak with one rotation i.e. if the invader
        # rotates into an upper strand, a leakage will not occur (on the upper strand). Variant of l_s_7 with longer sequences of domains.
        l_s_11 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L1 LA S* T^* RA R1}",
                                                        "{L' L2}<L LB>[S T^]<RB R>{R2 R'}")))[0].products.keys())
        exp_res_11 = {"{L' L2 S* T^* R2 R'}", "{L1 LA}<L LB>[S T^]<RB R>{RA R1}"}
        self.assertEqual(set(), set.difference(l_s_11, exp_res_11))

    def test_leakage_rule_yields_correctly_when_lower_strand_can_only_invade_as_upper_strand_long(self):
        # Test the LS rule when the invader strand is a lower strand which can only cause a leak if it rotates into an upper strand.
        # Variant of l_s_8 but with longer sequences of domains.
        l_s_12 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L1 LA T^ S RA R1}",
                                                        "{L' L2}<L LB>[S T^]<RB R>{R2 R'}")))[0].products.keys())
        exp_res_12 = {"<L LB S T^ RB R>", "{L' L2}<R1 RA>[S T^]<LA L1>{R2 R'}"}
        self.assertEqual(set(), set.difference(l_s_12, exp_res_12))

    def test_leakage_rule_does_not_displace_an_upper_strand_attached_to_a_previous_gate(self):
        # Test the LS rule does not displace an upper strand which connects directly to the previous gate.
        l_s_13 = set(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 S T R1>", "[A]<B>::{L'}<L>[S T]<R>{R'}"))))
        self.assertEqual(set(), l_s_13)

    def test_leakage_rule_does_not_displace_an_upper_strand_attached_to_a_following_gate(self):
        # Test the LS rule does not displace an upper strand which connects directly to the following gate.
        l_s_14 = set(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 S T R1>", "[A]<B>:{L'}<L>[S T]<R>{R'}::<C>[D]"))))
        self.assertEqual(set(), l_s_14)

    def test_leakage_rule_does_not_displace_a_lower_strand_attached_to_a_previous_gate(self):
        # Test the LS rule does not displace a lower strand which connects directly to the previous gate.
        l_s_15 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L1 S R1}", "[A]<B>:{L'}<L>[S^]<R>{R'}"))))
        self.assertEqual(set(), l_s_15)

    def test_leakage_rule_does_not_displace_a_lower_strand_attached_to_a_following_gate(self):
        # Test the LS rule does not displace an lower strand which connects directly to the following gate.
        l_s_16 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L1 S T R1}", "[A]<B>::{L'}<L>[S T]<R>{R'}:<C>[D]"))))
        self.assertEqual(set(), l_s_16)


class TestToeholdLeakageRule(unittest.TestCase):
    from stocal.examples.dsd import ToeholdLeakageRule
    Rule = ToeholdLeakageRule

    def test_lakin_l_t_example(self):
        # Test that the basic LT example from the Lakin paper can be replicated with the Leakage Rule.
        l_t_1 = set(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 S R1>", "{L'}<L>[S N^]<R>{R'}")))[0].products.keys())
        exp_res_1 = {"<L S N^ R>", "{L'}<L1>[S]<R1>{N^* R'}"}
        self.assertEqual(set(), set.difference(l_t_1, exp_res_1))

    def test_extended_lakin_l_t_example(self):
        # Test a different version of the LT example from the Lakin paper, with more domains on the double strand.
        l_t_2 = set(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 S K^ R1>", "{L'}<L>[S K^ N^]<R>{R'}")))[0].products.keys())
        exp_res_2 = {"<L S K^ N^ R>", "{L'}<L1>[S K^]<R1>{N^* R'}"}
        self.assertEqual(set(), set.difference(l_t_2, exp_res_2))

    def test_lower_strand_version_of_lakin_l_t_example(self):
        # Test that the basic (rotated) LT example from the Lakin paper can be replicated with the Leakage Rule.
        l_t_3 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L1 S* R1}", "{L'}<L>[S N^]<R>{R'}")))[0].products.keys())
        exp_res_3 = {"{L' S* N^* R'}", "{L1}<L>[S]<N^ R>{R1}"}
        self.assertEqual(set(), set.difference(l_t_3, exp_res_3))

    def test_extended_lower_strand_version_of_lakin_l_t_example(self):
        # Test that the basic (rotated) LT example from the Lakin paper can be replicated with the Leakage Rule.
        l_t_4 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L1 S* B^* R1}", "{L'}<L>[S B^ N^]<R>{R'}")))[0].products.keys())
        exp_res_4 = {"{L' S* B^* N^* R'}", "{L1}<L>[S B^]<N^ R>{R1}"}
        self.assertEqual(set(), set.difference(l_t_4, exp_res_4))

    def test_toehold_leak_where_upper_strand_only_initiates_leak_after_rotating_into_a_lower_strand(self):
        # Test that the basic LT example from the Lakin paper can be replicated, even when the strand is passed at the wrong rotation
        # and cannot initiate the leak until it rotates back to its original position.
        l_t_5 = set(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 S* R1>", "{L'}<L>[S N^]<R>{R'}")))[0].products.keys())
        exp_res_5 = {"{L' S* N^* R'}", "{R1}<L>[S]<N^ R>{L1}"}
        self.assertEqual(set(), set.difference(l_t_5, exp_res_5))

    def test_toehold_leak_where_lower_strand_only_initiates_leak_after_rotating_into_an_upper_strand(self):
        # Test that the basic LT example from the Lakin paper can be replicated, even when the strand is passed at the wrong rotation
        # and cannot initiate the leak until it rotates back to its original position.
        l_t_6 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{R1 S L1}", "{L'}<L>[S N^]<R>{R'}")))[0].products.keys())
        exp_res_6 = {"<L S N^ R>", "{L'}<L1>[S]<R1>{N^* R'}"}
        self.assertEqual(set(), set.difference(l_t_6, exp_res_6))

    def test_toehold_leak_with_toehold_at_start_of_double_strand_with_upper_invader_strand(self):
        # Test that the basic LT example from the Lakin paper can be replicated in reverse, right to left, when the
        # toehold occurs at the start of the double strand.
        l_t_7 = set(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 S R1>", "{L'}<L>[N^ S]<R>{R'}")))[0].products.keys())
        exp_res_7 = {"<L N^ S R>", "{L' N^*}<L1>[S]<R1>{R'}"}
        self.assertEqual(set(), set.difference(l_t_7, exp_res_7))

    def test_toehold_leak_with_toehold_at_start_of_double_strand_with_lower_invader_strand(self):
        # Test that the basic LT example from the Lakin paper can be replicated in reverse, right to left, when the
        # toehold occurs at the start of the double strand and the invader is a lower strand.
        l_t_8 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L1 S* R1}", "{L'}<L>[N^ S]<R>{R'}")))[0].products.keys())
        exp_res_8 = {"{L' N^* S* R'}", "{L1}<L N^>[S]<R>{R1}"}
        self.assertEqual(set(), set.difference(l_t_8, exp_res_8))

    def test_extended_lakin_l_t_example_with_toehold_at_start(self):
        # Test that the basic LT example from the Lakin paper can be replicated with the Leakage Rule.
        l_t_9 = set(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 N^ S R1>", "{L'}<L>[N^ S]<R>{R'}"))))
        exp_res_9 = {"<L N^ S R>", "{L' N^*}<L1>[S]<R1>{R'}"}
        self.assertEqual(set(), l_t_9)

    def test_lakin_l_s_example_does_not_yield_any_results_from_the_l_t_rule(self):
        # Test that the LT rule is not applied to the basic LS example from the Lakin paper.
        l_t_1 = set(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 S R1>", "{L'}<L>[S]<R>{R'}"))))
        self.assertEqual(set(), set.difference(l_t_1, set()))

    def test_that_a_rotated_lakin_l_s_example_does_not_yield_any_results_from_the_l_t_rule(self):
        # Test that the LT rule is not applied to the rotated (lower strand version) of the LS example from the Lakin paper.
        l_t_2 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L1 S R1}", "{L'}<L>[S]<R>{R'}"))))
        self.assertEqual(set(), set.difference(l_t_2, set()))

    def test_that_the_l_t_rule_does_not_apply_to_short_double_toeholds(self):
        # Test that the leakage rule does not yield any results when the short double strand has form [N^].
        l_t_3 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L1 S R1}", "{L'}<L>[S^]<R>{R'}"))))
        self.assertEqual(set(), l_t_3)

    def test_that_invader_strand_cannot_yield_a_toehold_leak_when_the_sequences_do_not_match(self):
        # Test that when the invader sequence of domains does not match the sequence of domains within the d_s of the
        # other input, no leakages are yielded
        l_t_4 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L1 A B^ C^ R1}", "{L'}<L>[A B C^]<R>{R'}"))))
        self.assertEqual(set(), l_t_4)

if __name__ == '__main__':
    unittest.main()
