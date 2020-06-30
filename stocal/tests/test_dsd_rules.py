"""Unit testing for rules in dsd.py """
import unittest
from stocal.tests.test_transitions import TestReactionRule as TestTransitionRule, TestMassAction


class TestBindingRule(unittest.TestCase):
    from stocal.examples.dsd import BindingRule
    Rule = BindingRule

    def test_lakin_rb_example_binds_correctly(self):
        # m_a_1 checks that the basic RB example from the Lakin paper yields the correct result.
        m_a_1 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L' N^* R'}", "<L N^ R>")))[0].products.keys())[0]
        self.assertEqual(m_a_1, "{L'}<L>[N^]<R>{R'}")

    def test_lakin_rb_example_binds_correctly_when_order_is_changed(self):
        # m_a_2 checks that m_a_1 yields the correct result even if the order is changed.
        m_a_2 = list(list(set(self.Rule.novel_reactions(self.Rule(), "<L N^ R>", "{L' N^* R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_2, "{L'}<L>[N^]<R>{R'}")

    def test_strand_to_strand_binding_when_binding_can_occur_in_multiple_places(self):
        # m_a_3 checks that the Binding Rule can yield multiple different bindings from the same inputs (if appropriate)
        m_a_3 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{S' N^* L' R'}", "<L N^ M N^>")))[0].products.keys())
        exp_res_3 = {"{S'}<L N^ M>[N^]{L' R'}", "{S'}<L>[N^]<M N^>{L' R'}"}
        self.assertEqual(set(), set.difference(m_a_3, exp_res_3))

    def test_strand_to_strand_binding_where_resulting_gate_has_no_lower_strand_preceding_the_double_strand(self):
        # Check slight variant of the Binding Rule, where the yielded result doesn't have a lower strand preceding the d_s
        m_a_4 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{ N^* L' R'}", "<L N^ M>")))[0].products.keys())[0]
        self.assertEqual(m_a_4, "<L>[N^]<M>{L' R'}")

    def test_strand_to_strand_binding_where_resulting_gate_has_no_lower_strand_after_the_double_strand(self):
        # Check slight variant of the Binding Rule, where the yielded result doesn't have a lower strand after the d_s
        m_a_5 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L' N^*}", "<L N^ M>")))[0].products.keys())[0]
        self.assertEqual(m_a_5, "{L'}<L>[N^]<M>")

    def test_strand_to_strand_binding_where_resulting_gate_has_no_upper_strand_preceding_the_double_strand(self):
        # Check slight variant of the Binding Rule, where the yielded result doesn't have an upper strand preceding the d_s
        m_a_6 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{A N^* L' R'}", "<N^ M>")))[0].products.keys())[0]
        self.assertEqual(m_a_6, "{A}[N^]<M>{L' R'}")

    def test_strand_to_strand_binding_where_resulting_gate_has_no_upper_strand_after_the_double_strand(self):
        # Check slight variant of the Binding Rule, where the yielded result doesn't have an upper strand after the d_s
        m_a_7 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L' N^* R}", "<L N^>")))[0].products.keys())[0]
        self.assertEqual(m_a_7, "{L'}<L>[N^]{R}")

    def test_simplest_binding_case(self):
        # Test the simplest strand to strand binding case, where the yielded result has just a single double toehold.
        m_a_8 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{N^*}", "<N^>")))[0].products.keys())[0]
        self.assertEqual(m_a_8, "[N^]")

    def test_lakin_fig_4_example_binds_correctly_strand_to_gate(self):
        # Check an example from Figure 4 of the Lakin paper
        m_a_9 = list(list(set(self.Rule.novel_reactions(self.Rule(), "<t^ x y>", "{t^*}[x]:[y u^]")))[0].products.keys())[0]
        self.assertEqual(m_a_9, "[t^]<x y>:[x]:[y u^]")

    def test_lakin_rp_example_binds_correctly_strand_to_gate(self):
        # m_a_1 checks that the basic RP example from the Lakin paper yields the correct result.
        m_a_10 = list(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 N^ S R1>", "{L' N^*}<L>[S R2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_10, "{L'}<L1>[N^]<S R1>:<L>[S R2]<R>{R'}")

    def test_binding_gate_to_gate_yields_no_results(self):
        # Check that binding does not occur between two gates.
        m_a_11 = set(self.Rule.novel_reactions(self.Rule(), "{N^* S' N^*}[C^]", "{L'}<L>[N^]<R>[M^]<S'>[A^]{B}"))
        self.assertEqual(m_a_11, set())

    def test_lower_strand_binding_to_gate(self):
        # Check that binding can occur between a lower strand and a gate.
        m_a_12 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{A C^*}", "{F}<B C^ G>[H^]<I>{J}")))[0].products.keys())[0]
        self.assertEqual(m_a_12, "{A}<B>[C^]::{F}<G>[H^]<I>{J}")

    def test_lower_strand_binding_to_gate_longer_sequence(self):
        # Check that binding can occur between a lower strand and a gate, and the gate is preceded by another gate.
        m_a_13 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{F}<B C^ D G>[H^]:{J K}<I L>[M^]<N>{O}", "{A C^* E}")))[0].products.keys())[0]
        self.assertEqual(m_a_13, "{A}<B>[C^]{E}::{F}<D G>[H^]:{J K}<I L>[M^]<N>{O}")

    def test_upper_strand_binding_to_gate(self):
        # Check that binding can occur between a lower strand and a gate.
        m_a_14 = list(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 N^ S R1>", "{L' N^*}<L>[S R2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_14, "{L'}<L1>[N^]<S R1>:<L>[S R2]<R>{R'}")


class TestUnbindingRule(TestTransitionRule):
    from stocal.examples.dsd import UnbindingRule
    Rule = UnbindingRule

    def test_lakin_ru_example_unbinds_correctly(self):
        # m_a_1 checks that the basic RU example from the Lakin paper yields the correct result.
        m_a_1 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[N^]<R>{R'}")))[0].products.keys())
        exp_res_1 = {"{L' N^* R'}", "<L N^ R>"}
        self.assertEqual(set(), set.difference(m_a_1, exp_res_1))

    def test_that_unbinding_generates_correct_results_on_longer_strands(self):
        # Check that RU correctly unbinds a gate which has more domains on its strands.
        m_a_2 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{B}<A>[D^]<C^ F>{C^* G}")))[0].products.keys())
        exp_res_2 = {"<A D^ C^ F>", "{B D^* C^* G}"}
        self.assertEqual(set(), set.difference(m_a_2, exp_res_2))

    def test_that_the_second_gate_in_a_system_can_unbind_correctly(self):
        # Test a system which consists of two gates, with one possible point of unbinding, on the second gate.
        m_a_3 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L1>[N^]<S R1>:<L>[S R2]<R>{R'}")))[0].products.keys())
        exp_res_3 = {"<L1 N^ S R1>", "{L' N^*}<L>[S R2]<R>{R'}"}
        self.assertEqual(set(), set.difference(m_a_3, exp_res_3))

    def test_a_system_which_can_unbind_correctly_at_multiple_points(self):
        # Test a system which can unbind at 3 different points.
        m_a_4 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{A}<B>[C^]<D>{E}::{F}<G>[H^]<I>{J}::{K}<L>[M^]<N>{O}")))[0].products.keys())
        exp_res_4 = {"{F}<B C^ D G>[H^]{J}::{K}<I L>[M^]<N>{O}", "{A C^* E}", "{A}<B>[C^]{E}::{K}<D G H^ I L>[M^]<N>{O}", "{F H^* J}",
                     "{A}<B>[C^]{E}::{F}<D G>[H^]<I L M^ N>{J}", "{K M^* O}"}
        self.assertEqual(set(), set.difference(m_a_4, exp_res_4))


class TestCoveringRule(TestTransitionRule):
    from stocal.examples.dsd import CoveringRule
    Rule = CoveringRule
    def test_lakin_rc_example_generates_correct_result(self):
        # m_a_1 checks that the basic RC  example from the Lakin paper yields the correct result.
        m_a_1 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S]<N^ R>{N^* R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_1, "{L'}<L>[S N^]<R>{R'}")

    def test_lakin_rc_example_generates_correct_result_right_to_left(self):
        # m_a_2 checks that the RC example works right to left, as well as left to right.
        m_a_2 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L' N^*}<L N^>[S]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_2, "{L'}<L>[N^ S]<R>{R'}")

    def test_covering_rule_left_to_right(self):
        # Test a basic variant of the covering rule RC, applied left to right.
        m_a_3 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[S]<N^ R>{N^* R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_3, "[S N^]<R>{R'}")

    def test_covering_rule_right_to_left(self):
        # Test a basic variant of the covering rule RC, applied right to left.
        m_a_4 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{R' N^*}<R N^>[S]")))[0].products.keys())[0]
        self.assertEqual(m_a_4, "{R'}<R>[N^ S]")

    def test_covering_rule_left_to_right_across_gates_joined_at_upper_strand(self):
        # Test the application of the covering rule across gates, left to right, where the gates are joined by an upper strand.
        m_a_5 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{A}<B>[C]{E^*}::{F}<E^ D>[G]")))[0].products.keys())[0]
        self.assertEqual(m_a_5, "{A}<B>[C E^]::{F}<D>[G]")
        # N.B: No right_to_left version of this exists, due to the chosen normal form.

    def test_covering_rule_left_to_right_across_gates_joined_at_upper_strand_variant(self):
        # A variation of the last test, where the lower domain which is being bound to is followed by other domains.
        m_a_6 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{A}<B>[C]{E^* Z}::{F}<E^ D>[G]")))[0].products.keys())[0]
        self.assertEqual(m_a_6, "{A}<B>[C E^]{Z}::{F}<D>[G]")
        # N.B: No right_to_left version of this exists, due to the chosen normal form.

    def test_covering_rule_left_to_right_variant(self):
        # Tests a variation of the covering rule where the gate which is being 'covered' is followed by another gate.
        m_a_7 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S]<N^ R>{N^* R'}::[A B]")))[0].products.keys())[0]
        self.assertEqual(m_a_7, "{L'}<L>[S N^]<R>{R'}::[A B]")
        # N.B: No right_to_left version of this exists, due to the chosen normal form.

    def test_covering_rule_left_to_right_variant_2(self):
        # Tests a variation of the covering rule where the gate which is being 'covered' is in between other gates.
        m_a_8 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[C D]<A>:{L'}<L>[S]<N^ R>{N^* R'}::[A B]")))[0].products.keys())[0]
        self.assertEqual(m_a_8, "[C D]<A>:{L'}<L>[S N^]<R>{R'}::[A B]")
        # N.B: No right_to_left version of this exists, due to the chosen normal form.



class TestMigrationRule(TestTransitionRule):
    from stocal.examples.dsd import MigrationRule
    Rule = MigrationRule

    def test_lakin_rm_example_upper_l_to_r(self):
        # m_a_1 checks that the basic RM example from the Lakin paper yields the correct result.
        m_a_1 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S R2>:<L1>[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_1, "{L'}<L>[S1 S]<R2>:<L1 S>[S2]<R>{R'}")

    def test_lakin_rm_example_lower_l_to_r(self):
        # Test variants of m_a_1 but when the overhang is on the lower strand:
        m_a_2 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S R2}::{L1}[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_2, "{L'}<L>[S1 S]{R2}::{L1 S}[S2]<R>{R'}")

    def test_lakin_rm_example_upper_l_to_r_second_overhang_only_in_result(self):
        # Check variant of m_a_1 where R2 is missing (so the result only has one overhang):
        m_a_3 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S>:<L1>[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_3, "{L'}<L>[S1 S]:<L1 S>[S2]<R>{R'}")

    def test_lakin_rm_example_lower_l_to_r_second_overhang_only_in_result(self):
        # Check lower_strand variant of m_a_1 where R2 is missing (so the result only has one overhang):
        m_a_4 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S}::{L1}[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_4, "{L'}<L>[S1 S]::{L1 S}[S2]<R>{R'}")

    def test_lakin_rm_example_upper_l_to_r_input_only_has_first_overhang(self):
        # Check variant of m_a_1 where the input only has the first overhang (i.e. L1 is missing)
        m_a_5 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S R2>:[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_5, "{L'}<L>[S1 S]<R2>:<S>[S2]<R>{R'}")

    def test_lakin_rm_example_lower_l_to_r_input_only_has_first_overhang(self):
        # Check lower strand variant of m_a_1 where the input only has the first overhang (i.e. L1 is missing)
        m_a_6 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S R2}::[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_6, "{L'}<L>[S1 S]{R2}::{S}[S2]<R>{R'}")

    def test_lakin_rm_example_upper_l_to_r_input_only_has_the_first_overhang_and_result_only_has_second_overhang(self):
        # Check variants of m_a_1 where R2 and L1 are missing:
        m_a_7 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S>:[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_7, "{L'}<L>[S1 S]:<S>[S2]<R>{R'}")

    def test_lakin_rm_example_lower_l_to_r_input_only_has_the_first_overhang_and_result_only_has_second_overhang(self):
        # Check variants of m_a_1 where R2 and L1 are missing:
        m_a_8 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S}::[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_8, "{L'}<L>[S1 S]::{S}[S2]<R>{R'}")

    def test_that_migration_rule_is_not_applied_to_lakin_displacement_example_rd(self):
        # Check that RM is not applied on the RD example, as the two should be mutually exclusive.
        m_a_9 = set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S R>:<L2>[S]<R2>{R'}"))
        self.assertEqual(m_a_9, set())

    def test_that_migration_rule_is_not_applied_to_lower_strand_version_of_lakin_displacement_example_rd(self):
        # Check that RM is not applied to the lower strand version of the RD example, as the rules should be mutually exclusive.
        m_a_10 = set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S R}::{L2}[S]<R2>{R'}"))
        self.assertEqual(m_a_10, set())

    def test_that_migration_rule_is_not_applied_to_lakin_displacement_example_fig_4a(self):
        # Check that the RM rule is not applied to the RD example from Figure 4a).
        m_a_11 = set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x]:[y u^]"))
        self.assertEqual(m_a_11, set())

    def test_that_migration_rule_is_not_applied_to_lower_strand_version_of_lakin_displacement_example_fig_4a(self):
        # Check that the RM rule is not applied to the lower strand version of the RD example from Figure 4a).
        m_a_12 = set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::[x]::[y u^]"))
        self.assertEqual(m_a_12, set())

    def test_upper_l_to_r_lakin_fig_4a_migration_example_correct(self):
        # Check the migration rule is applied correctly to the example from Figure 4a) of Lakin's paper.
        m_a_13 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]<y>:[y u^]")))[0].products.keys())[0]
        self.assertEqual(m_a_13, "[t^ x y]:<y>[u^]")

    def test_lower_l_to_r_lakin_fig_4a_migration_example_correct(self):
        # Check the migration rule is applied correctly to the lower strand version of the example from Figure 4a) of Lakin's paper.
        m_a_14 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]{y}::[y u^]")))[0].products.keys())[0]
        self.assertEqual(m_a_14, "[t^ x y]::{y}[u^]")

    def test_migration_rule_upper_l_to_r_variant_1(self):
        # Check example where the second gate involved in migration is connected to a third gate via the upper strand.
        m_a_15 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x v]::[y u^]")))[0].products.keys())[0]
        self.assertEqual("[t^ x]<y>:<x>[v]::[y u^]", m_a_15)

    def test_migration_rule_lower_l_to_r_variant_1(self):
        # Check example where the second gate involved in migration is connected to a third gate via the lower strand.
        m_a_16 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::[x v]:[y u^]")))[0].products.keys())[0]
        self.assertEqual("[t^ x]{y}::{x}[v]:[y u^]", m_a_16)

    def test_lakin_rm_example_upper_r_to_l(self):
        # m_a_17 checks that the basic RM example from the Lakin paper yields the correct result - when done in reverse (right to left).
        m_a_17 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1 S]<R2>:<L1 S>[S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_17, "{L'}<L>[S1]<S R2>:<L1>[S S2]<R>{R'}")

    def test_lakin_rm_example_lower_r_to_l(self):
        # m_a_17 checks that the lower strand version of the RM example from the Lakin paper can be performed left to right (reverse)
        m_a_18 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1 S]<R2>:<L1 S>[S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_18, "{L'}<L>[S1]<S R2>:<L1>[S S2]<R>{R'}")

    def test_lakin_rm_example_upper_r_to_l_second_overhang_only_in_input(self):
        # Check variant of RM (applied right to left) where the input only has the second overhang
        m_a_19 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1 S]:<L1 S>[S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_19, "{L'}<L>[S1]<S>:<L1>[S S2]<R>{R'}")

    def test_lakin_rm_example_lower_r_to_l_second_overhang_only_in_input(self):
        # Check lower strand variant of RM (applied right to left) where the input only has the second overhang
        m_a_20 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1 S]::{L1 S}[S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_20, "{L'}<L>[S1]{S}::{L1}[S S2]<R>{R'}")

    def test_lakin_rm_example_upper_r_to_l_result_only_has_first_overhang(self):
        # Check variant of RM (applied right to left) where the result only has the first overhang
        m_a_21 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1 S]<R2>:<S>[S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_21, "{L'}<L>[S1]<S R2>:[S S2]<R>{R'}")

    def test_lakin_rm_example_lower_r_to_l_result_only_has_first_overhang(self):
        # Check lower strand variant of RM (appied right to left) where the result only has the first overhang
        m_a_22 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1 S]{R2}::{S}[S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_22, "{L'}<L>[S1]{S R2}::[S S2]<R>{R'}")

    def test_lakin_rm_example_upper_r_to_l_input_only_has_the_second_overhang_and_result_only_has_first_overhang(self):
        # Check variant of Lakin's RM rule (applied right to left) where R2 and L1 are missing:
        m_a_23 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1 S]:<S>[S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_23, "{L'}<L>[S1]<S>:[S S2]<R>{R'}")

    def test_lakin_rm_example_lower_r_to_l_input_only_has_the_second_overhang_and_result_only_has_first_overhang(self):
        # Check lower strand variant of Lakin's RM rule (applied right to left) where R2 and L1 are missing:
        m_a_24 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1 S]::{S}[S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_24, "{L'}<L>[S1]{S}::[S S2]<R>{R'}")

    def test_upper_r_to_r_lakin_fig_4a_migration_example_correct(self):
        # Check the migration rule is applied correctly (in reverse) to the example from Figure 4a) of Lakin's paper.
        m_a_25 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x y]:<y>[u^]")))[0].products.keys())[0]
        self.assertEqual(m_a_25, "[t^ x]<y>:[y u^]")

    def test_lower_r_to_l_lakin_fig_4a_migration_example_correct(self):
        # Check the migration rule is applied correctly (in reverse) to the lower strand version of the Figure 4a eg. in Lakin's paper.
        m_a_26 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x y]::{y}[u^]")))[0].products.keys())[0]
        self.assertEqual(m_a_26, "[t^ x]{y}::[y u^]")

    def test_migration_rule_upper_r_to_l_variant_1(self):
        # Check right to left example where the second gate involved in migration is connected to a third gate via the upper strand.
        m_a_27 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]<y>:<x>[v]::[y u^]")))[0].products.keys())[0]
        self.assertEqual("[t^]<x y>:[x v]::[y u^]", m_a_27)

    def test_migration_rule_lower_r_to_l_variant_1(self):
        # Check right to left example where the second gate involved in migration is connected to a third gate via the lower strand.
        m_a_28 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]{y}::{x}[v]:[y u^]")))[0].products.keys())[0]
        self.assertEqual("[t^]{x y}::[x v]:[y u^]", m_a_28)


        #
        # # Test variants of m_a_1 but when the overhang is on the lower strand:
        # m_a_4 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S R2}::{L1}[S S2]<R>{R'}")))[0].products.keys())[0]
        # self.assertEqual(m_a_4, "{L'}<L>[S1 S]{R2}::{L1 S}[S2]<R>{R'}")
        #
        # # Check lower strand equivalents of m_a_2 and m_a_3
        # m_a_5 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S}::{L1}[S S2]<R>{R'}")))[0].products.keys())[0]
        # self.assertEqual(m_a_5, "{L'}<L>[S1 S]::{L1 S}[S2]<R>{R'}")
        # m_a_6 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S}::[S S2]<R>{R'}")))[0].products.keys())[0]
        # self.assertEqual(m_a_6, "{L'}<L>[S1 S]::{S}[S2]<R>{R'}")
        #
        # # Check that RM is not applied on the RD example, as the two should be mutually exclusive.
        # m_a_7 = set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S R>:<L2>[S]<R2>{R'}"))
        # self.assertEqual(m_a_7, set())
        # m_a_8 = set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S R}::{L2}[S]<R2>{R'}"))
        # self.assertEqual(m_a_8, set())
        #
        # # Check that the RM rule is not applied to the RD example from Figure 4a).
        # m_a_9 = set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x]:[y u^]"))
        # self.assertEqual(m_a_9, set())
        # m_a_10 = set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::[x]::[y u^]"))
        # self.assertEqual(m_a_10, set())
        #
        # # Check the migration rule is applied correctly to the example from Figure 4a) of Lakin's paper.
        # m_a_11 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]<y>:[y u^]")))[0].products.keys())[0]
        # self.assertEqual(m_a_11, "[t^ x y]:<y>[u^]")
        # m_a_12 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]{y}::[y u^]")))[0].products.keys())[0]
        # self.assertEqual(m_a_12, "[t^ x y]::{y}[u^]")
        #
        # # Check additional variants:
        # m_a_13 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x v]::[y u^]")))[0].products.keys())[0]
        # self.assertEqual("[t^ x]<y>:<x>[v]::[y u^]", m_a_13)


class TestDisplacementRule(TestTransitionRule):
    from stocal.examples.dsd import DisplacementRule
    Rule = DisplacementRule

    def test_reduction_rule_generates_a_b_c(self):
        # Test the rule reduction example RD from Lakin's paper.
        m_a_1 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S R>:<L2>[S]<R2>{R'}")))[0].products.keys())
        exp_res_1 = {"<L2 S R2>", "{L'}<L>[S1 S]<R>{R'}"}
        self.assertEqual(set(), set.difference(m_a_1, exp_res_1))

        # Test the lower strand equivalent of the reduction example RD from Lakin's paper.
        m_a_2 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S R}::{L2}[S]<R2>{R'}")))[0].products.keys())
        exp_res_2 = {"{L2 S R'}", "{L'}<L>[S1 S]<R2>{R}"}
        self.assertEqual(set(), set.difference(m_a_2, exp_res_2))

        # m_a_3 and m_a_4 checks that the application of the Reduction rule from Figure 4 works as expected.
        m_a_3 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x]:[y u^]")))[0].products.keys())
        exp_res_3 = {"<x>", "[t^ x]<y>:[y u^]"}
        self.assertEqual(set(), set.difference(m_a_3, exp_res_3))
        m_a_4 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::[x]::[y u^]")))[0].products.keys())
        exp_res_4 = {"{x}", "[t^ x]{y}::[y u^]"}
        self.assertEqual(set(), set.difference(m_a_4, exp_res_4))

        # m_a_5 checks that the Displacement rule does not get applied to the Migration example from Figure 4a of the Lakin paper.
        m_a_5 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]<y>:[y u^]"))))
        self.assertEqual(set(), m_a_5)
        # m_a_6 checks that the example from 4a with flipped orientation cannot yield displacement products.
        m_a_6 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]{y}::[y u^]"))))
        self.assertEqual(m_a_6, set())

        # Check that other systems where migration can occur cannot be displaced:
        m_a_7 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x v]::[y u^]"))))
        self.assertEqual(set(), m_a_7)
        m_a_8 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::[x v]:[y u^]"))))
        self.assertEqual(set(), m_a_8)

        # This test checks that applying the displacement rule along an upper strand works, when the toehold which is being
        # displaced is connected along its lower strand to the next gate (left to right).
        m_a_9 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x]::[y u^]")))[0].products.keys())
        exp_res_9 = {"[t^ x]<y>", "<x>[y u^]"}
        self.assertEqual(set(), set.difference(m_a_9, exp_res_9))
        # This is a variant of m_a_6 but switching orientation.
        m_a_10 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::[x]:[y u^]")))[0].products.keys())
        exp_res_10 = {"[t^ x]{y}", "{x}[y u^]"}
        self.assertEqual(set(), set.difference(m_a_10, exp_res_10))

        # More variants of m_a_9 and m_a_10:
        m_a_11 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:<R>[x]::[y u^]")))[0].products.keys())
        exp_res_11 = {"[t^ x]<y>", "<R x>[y u^]"}
        self.assertEqual(set(), set.difference(m_a_11, exp_res_11))
        m_a_12 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::{R}[x]:[y u^]")))[0].products.keys())
        exp_res_12 = {"[t^ x]{y}", "{R x}[y u^]"}
        self.assertEqual(set(), set.difference(m_a_12, exp_res_12))

        m_a_13 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:<r>[x]{g}::[y u^]")))[0].products.keys())
        exp_res_13 = {"[t^ x]<y>{g}", "<r x>[y u^]"}
        self.assertEqual(set(), set.difference(m_a_13, exp_res_13))
        m_a_14 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::{r}[x]<g>:[y u^]")))[0].products.keys())
        exp_res_14 = {"[t^ x]<g>{y}", "{r x}[y u^]"}
        self.assertEqual(set(), set.difference(m_a_14, exp_res_14))


if __name__ == '__main__':
    unittest.main()
