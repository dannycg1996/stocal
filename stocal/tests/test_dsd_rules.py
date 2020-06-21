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

    def test_strand_to_strand_binding_generates_a_b_c(self):

        #Test the simplest unbinding cases, where the system consists of a single gate.
        self.test_single_unbinding("{L'}<L>[N^]<R>{R'}", "{L' N^* R'}", "<L N^ R>")
        self.test_single_unbinding("{B}<A>[D^]<C^ F>{C^* G}", "<A D^ C^ F>", "{B D^* C^* G}")

    def test_single_unbinding(self, whole, part_1, part_2):
        """This function takes a system (which should only contain one double toehold which can unbind) and two expected
        results (from being unbound) and then checks that the rule works"""
        result_1 = list(list(set(self.Rule.novel_reactions(self.Rule(), whole)))[0].products.keys())[0]
        result_2 = list(list(set(self.Rule.novel_reactions(self.Rule(), whole)))[0].products.keys())[1]
        if result_1 == part_1:
            self.assertEqual(result_1, part_1)
            self.assertEqual(result_2, part_2)
        else:
            self.assertEqual(result_1, part_2)
            self.assertEqual(result_2, part_1)


        # TODO: Can I unbind [A^ B^]?





if __name__ == '__main__':
    unittest.main()
