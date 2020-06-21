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

    # TODO: Can I unbind [A^ B^]?
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


if __name__ == '__main__':
    unittest.main()
