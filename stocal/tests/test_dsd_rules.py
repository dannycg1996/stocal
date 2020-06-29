"""Bugs reported in github issue tracker
"""
import unittest
from stocal.tests.test_transitions import TestReactionRule as TestTransitionRule, TestMassAction


class TestBindingRule(TestTransitionRule):
    from stocal.examples.dsd import BindingRule
    Rule = BindingRule

    def test_strand_to_strand_binding_generates_a_b_c(self):
        # m_a_1 checks that the basic RB example from the Lakin paper yields the correct result.
        m_a_1 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L' N^* R'}", "<L N^ R>")))[0].products.keys())[0]
        self.assertEqual(m_a_1, "{L'}<L>[N^]<R>{R'}")

        # m_a_2 checks that m_a_1 yields the correct result even if the order is changed.
        m_a_2 = list(list(set(self.Rule.novel_reactions(self.Rule(), "<L N^ R>", "{L' N^* R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_2, "{L'}<L>[N^]<R>{R'}")

        # m_a_3 and m_a_4 check that the Binding Rule can yield multiple different bindings from the same inputs (if appropriate)
        m_a_3 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{S' N^* L' R'}", "<L N^ M N^>")))[0].products.keys())[0]
        m_a_4 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{S' N^* L' R'}", "<L N^ M N^>")))[1].products.keys())[0]
        if m_a_3 == "{S'}<L N^ M>[N^]{L' R'}":
            self.assertEqual(m_a_3, "{S'}<L N^ M>[N^]{L' R'}")
            self.assertEqual(m_a_4, "{S'}<L>[N^]<M N^>{L' R'}")
        else:
            self.assertEqual(m_a_3, "{S'}<L>[N^]<M N^>{L' R'}")
            self.assertEqual(m_a_4, "{S'}<L N^ M>[N^]{L' R'}")

        # Check slight variants of the Binding Rule, where the yielded result doesn't all 5 strands i.e. {}<>[]<>{}
        m_a_5 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{ N^* L' R'}", "<L N^ M>")))[0].products.keys())[0]
        self.assertEqual(m_a_5, "<L>[N^]<M>{L' R'}")
        m_a_6 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{N^*}", "<N^>")))[0].products.keys())[0]
        self.assertEqual(m_a_6, "[N^]")

        # Check an example from Figure 4 of the Lakin paper
        m_a_7 = list(list(set(self.Rule.novel_reactions(self.Rule(), "<t^ x y>", "{t^*}[x]:[y u^]")))[0].products.keys())[0]
        self.assertEqual(m_a_7, "[t^]<x y>:[x]:[y u^]")

    def test_strand_to_gate_binding_generates_a_b_c(self):
        # m_a_1 checks that the basic RP example from the Lakin paper yields the correct result.
        m_a_1 = list(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 N^ S R1>", "{L' N^*}<L>[S R2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_1, "{L'}<L1>[N^]<S R1>:<L>[S R2]<R>{R'}")

        # Check that binding does not occur between two gates.
        m_a_2 = set(self.Rule.novel_reactions(self.Rule(), "{N^* S' N^*}[C^]", "{L'}<L>[N^]<R>[M^]<S'>[A^]{B}"))
        self.assertEqual(m_a_2, set())

        # Strand to gate variants
        m_a_3 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{A C^*}", "{F}<B C^ G>[H^]<I>{J}")))[0].products.keys())[0]
        self.assertEqual(m_a_3, "{A}<B>[C^]::{F}<G>[H^]<I>{J}")
        m_a_4 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{F}<B C^ D G>[H^]:{J K}<I L>[M^]<N>{O}", "{A C^* E}")))[0].products.keys())[0]
        self.assertEqual(m_a_4, "{A}<B>[C^]{E}::{F}<D G>[H^]:{J K}<I L>[M^]<N>{O}")
        m_a_5 = list(list(set(self.Rule.novel_reactions(self.Rule(), "<L1 N^ S R1>", "{L' N^*}<L>[S R2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_5, "{L'}<L1>[N^]<S R1>:<L>[S R2]<R>{R'}")


class TestUnbindingRule(TestTransitionRule):
    from stocal.examples.dsd import UnbindingRule
    Rule = UnbindingRule

    def test_strand_to_strand_binding_generates_a_b_c(self):
        # Test the simplest unbinding cases, where the system consists of a single gate.
        m_a_1 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[N^]<R>{R'}")))[0].products.keys())
        exp_res_1 = {"{L' N^* R'}", "<L N^ R>"}
        self.assertEqual(set(), set.difference(m_a_1,exp_res_1))

        m_a_2 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{B}<A>[D^]<C^ F>{C^* G}")))[0].products.keys())
        exp_res_2 = {"<A D^ C^ F>", "{B D^* C^* G}"}
        self.assertEqual(set(), set.difference(m_a_2,exp_res_2))

        # Test a system which consists of two gates, with one possible point of unbinding.
        m_a_3 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L1>[N^]<S R1>:<L>[S R2]<R>{R'}")))[0].products.keys())
        exp_res_3 = {"<L1 N^ S R1>", "{L' N^*}<L>[S R2]<R>{R'}"}
        self.assertEqual(set(), set.difference(m_a_3,exp_res_3))

        # Test a system which can unbind at 3 different points.
        m_a_4 = set(list(set(self.Rule.novel_reactions(self.Rule(),"{A}<B>[C^]<D>{E}::{F}<G>[H^]<I>{J}::{K}<L>[M^]<N>{O}")))[0].products.keys())
        exp_res_4 = {"{F}<B C^ D G>[H^]{J}::{K}<I L>[M^]<N>{O}", "{A C^* E}", "{A}<B>[C^]{E}::{K}<D G H^ I L>[M^]<N>{O}","{F H^* J}",
                     "{A}<B>[C^]{E}::{F}<D G>[H^]<I L M^ N>{J}", "{K M^* O}"}
        self.assertEqual(set(), set.difference(m_a_4,exp_res_4))


class TestCoveringRule(TestTransitionRule):
    from stocal.examples.dsd import CoveringRule
    Rule = CoveringRule

    def test_covering_rule_generates_a_b_c(self):
        # m_a_1 checks that the basic RC  example from the Lakin paper yields the correct result.
        m_a_1 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S]<N^ R>{N^* R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_1, "{L'}<L>[S N^]<R>{R'}")
        # m_a_2 checks that the RC example works right to left, as well as left to right.
        m_a_2 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L' N^*}<L N^>[S]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_2, "{L'}<L>[N^ S]<R>{R'}")

        # Check variants:
        m_a_3 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[S]<N^ R>{N^* R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_3, "[S N^]<R>{R'}")
        m_a_4 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{A}<B>[C]{E^*}::{F}<E^ D>[G]")))[0].products.keys())[0]
        self.assertEqual(m_a_4, "{A}<B>[C E^]::{F}<D>[G]")
        m_a_5 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{A}<B>[C]{E^* Z}::{F}<E^ D>[G]")))[0].products.keys())[0]
        self.assertEqual(m_a_5, "{A}<B>[C E^]{Z}::{F}<D>[G]")
        m_a_6 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{A}<B>[C D]<E^ F>:{E^* G}[H]")))[0].products.keys())[0]
        self.assertEqual(m_a_6, "{A}<B>[C D E^]<F>:{G}[H]")
        m_a_7 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S]<N^ R>{N^* R'}::[A B]")))[0].products.keys())[0]
        self.assertEqual(m_a_7, "{L'}<L>[S N^]<R>{R'}::[A B]")
        m_a_8 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[C D]<A>:{L'}<L>[S]<N^ R>{N^* R'}::[A B]")))[0].products.keys())[0]
        self.assertEqual(m_a_8, "[C D]<A>:{L'}<L>[S N^]<R>{R'}::[A B]")


class TestMigrationRule(TestTransitionRule):
    from stocal.examples.dsd import MigrationRule
    Rule = MigrationRule

    def test_migration_rule_generates_a_b_c(self):
        # m_a_1 checks that the basic RM example from the Lakin paper yields the correct result.
        m_a_1 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S R2>:<L1>[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_1, "{L'}<L>[S1 S]<R2>:<L1 S>[S2]<R>{R'}")

        # Check variants of m_a_1 where R2 and L1 are missing:
        m_a_2 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S>:<L1>[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_2, "{L'}<L>[S1 S]:<L1 S>[S2]<R>{R'}")
        m_a_3 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S>:[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_3, "{L'}<L>[S1 S]:<S>[S2]<R>{R'}")

        # Test variants of m_a_1 but when the overhang is on the lower strand:
        m_a_4 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S R2}::{L1}[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_4, "{L'}<L>[S1 S]{R2}::{L1 S}[S2]<R>{R'}")

        # Check lower strand equivalents of m_a_2 and m_a_3
        m_a_5 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S}::{L1}[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_5, "{L'}<L>[S1 S]::{L1 S}[S2]<R>{R'}")
        m_a_6 = list(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S}::[S S2]<R>{R'}")))[0].products.keys())[0]
        self.assertEqual(m_a_6, "{L'}<L>[S1 S]::{S}[S2]<R>{R'}")

        # Check that RM is not applied on the RD example, as the two should be mutually exclusive.
        m_a_7 = set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S R>:<L2>[S]<R2>{R'}"))
        self.assertEqual(m_a_7, set())
        m_a_8 = set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S R}::{L2}[S]<R2>{R'}"))
        self.assertEqual(m_a_8, set())

        # Check that the RM rule is not applied to the RD example from Figure 4a).
        m_a_9 = set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x]:[y u^]"))
        self.assertEqual(m_a_9, set())
        m_a_10 = set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::[x]::[y u^]"))
        self.assertEqual(m_a_10, set())

        # Check the migration rule is applied correctly to the example from Figure 4a) of Lakin's paper.
        m_a_11 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]<y>:[y u^]")))[0].products.keys())[0]
        self.assertEqual(m_a_11, "[t^ x y]:<y>[u^]")
        m_a_12 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]{y}::[y u^]")))[0].products.keys())[0]
        self.assertEqual(m_a_12, "[t^ x y]::{y}[u^]")

        # Check additional variants:
        m_a_13 = list(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x v]::[y u^]")))[0].products.keys())[0]
        self.assertEqual("[t^ x]<y>:<x>[v]::[y u^]" ,m_a_13)


class TestDisplacementRule(TestTransitionRule):
    from stocal.examples.dsd import DisplacementRule
    Rule = DisplacementRule

    def test_reduction_rule_generates_a_b_c(self):
        # Test the rule reduction example RD from Lakin's paper.
        m_a_1 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]<S R>:<L2>[S]<R2>{R'}")))[0].products.keys())
        exp_res_1 = {"<L2 S R2>", "{L'}<L>[S1 S]<R>{R'}"}
        self.assertEqual(set(), set.difference(m_a_1,exp_res_1))

        # Test the lower strand equivalent of the reduction example RD from Lakin's paper.
        m_a_2 = set(list(set(self.Rule.novel_reactions(self.Rule(), "{L'}<L>[S1]{S R}::{L2}[S]<R2>{R'}")))[0].products.keys())
        exp_res_2 = {"{L2 S R'}", "{L'}<L>[S1 S]<R2>{R}"}
        self.assertEqual(set(), set.difference(m_a_2,exp_res_2))

        # m_a_3 and m_a_4 checks that the application of the Reduction rule from Figure 4 works as expected.
        m_a_3 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x]:[y u^]")))[0].products.keys())
        exp_res_3 = {"<x>", "[t^ x]<y>:[y u^]"}
        self.assertEqual(set(), set.difference(m_a_3,exp_res_3))
        m_a_4 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::[x]::[y u^]")))[0].products.keys())
        exp_res_4 = {"{x}", "[t^ x]{y}::[y u^]"}
        self.assertEqual(set(), set.difference(m_a_4,exp_res_4))

        # m_a_5 checks that the Displacement rule does not get applied to the Migration example from Figure 4a of the Lakin paper.
        m_a_5 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]<y>:[y u^]"))))
        self.assertEqual(set(), m_a_5)
        # m_a_6 checks that the example from 4a with flipped orientation cannot yield displacement products.
        m_a_6 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^ x]{y}::[y u^]"))))
        self.assertEqual(m_a_6, set())

        # Check that other systems where migration can occur cannot be displaced:
        m_a_7 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x v]::[y u^]"))))
        self.assertEqual(set(),m_a_7)
        m_a_8 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::[x v]:[y u^]"))))
        self.assertEqual(set(),m_a_8)

        # This test checks that applying the displacement rule along an upper strand works, when the toehold which is being
        # displaced is connected along its lower strand to the next gate (left to right).
        m_a_9 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:[x]::[y u^]")))[0].products.keys())
        exp_res_9 = {"[t^ x]<y>","<x>[y u^]"}
        self.assertEqual(set(), set.difference(m_a_9,exp_res_9))
        # This is a variant of m_a_6 but switching orientation.
        m_a_10 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::[x]:[y u^]")))[0].products.keys())
        exp_res_10 = {"[t^ x]{y}","{x}[y u^]"}
        self.assertEqual(set(), set.difference(m_a_10,exp_res_10))

        # More variants of m_a_9 and m_a_10:
        m_a_11 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:<R>[x]::[y u^]")))[0].products.keys())
        exp_res_11 = {"[t^ x]<y>", "<R x>[y u^]"}
        self.assertEqual(set(),set.difference(m_a_11,exp_res_11))
        m_a_12 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::{R}[x]:[y u^]")))[0].products.keys())
        exp_res_12 = {"[t^ x]{y}", "{R x}[y u^]"}
        self.assertEqual(set(),set.difference(m_a_12,exp_res_12))

        m_a_13 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]<x y>:<r>[x]{g}::[y u^]")))[0].products.keys())
        exp_res_13 = {"[t^ x]<y>{g}", "<r x>[y u^]"}
        self.assertEqual(set(),set.difference(m_a_13,exp_res_13))
        m_a_14 = set(list(set(self.Rule.novel_reactions(self.Rule(), "[t^]{x y}::{r}[x]<g>:[y u^]")))[0].products.keys())
        exp_res_14 = {"[t^ x]<g>{y}", "{r x}[y u^]"}
        self.assertEqual(set(),set.difference(m_a_14,exp_res_14))


if __name__ == '__main__':
    unittest.main()
