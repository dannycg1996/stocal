import unittest
from Main import regex_match

class MyTestCase(unittest.TestCase):
    def test_lower_th_regex(self):
        self.assertEqual(regex_match("{A' 1 B^ F G}", "lower_th"),['B'])
        self.assertEqual(regex_match("{L' A^* R'}{L' B^* R'} {L' C^ R'}", "lower_th"), ['C'])#
        self.assertEqual(regex_match("<A' 1 B^ F G>", "lower_th"), [])

    def test_upper_th_regex(self):
        self.assertEqual(regex_match("<L D^ R>", "upper_th"), ['D'])
        self.assertEqual(regex_match("{L D^ R}", "upper_th"), [])
        self.assertEqual(regex_match("{1 2 F^}[3 4 5]<L D^ R>{S U V}<N X^ M>", "upper_th"), ['D','X'])
        self.assertEqual(regex_match("{1 2 F^}[3 4 5]<L D^ R>{S U V}<N X^* M>", "upper_th"), ['D'])

    def test_lower_thc_regex(self):
        self.assertEqual(regex_match("{L' A^* }", "lower_th_c"),['A'])
        self.assertEqual(regex_match("{L' A^* R'}{L' B^* R'} {L' C^ R'}", "lower_th_c"),['A','B'])
        self.assertEqual(regex_match("{L' A^* R'}<1 2^ 3>{L' C^ R'}", "lower_th_c"),['A'])
        self.assertEqual(regex_match("<1 2^ 3>{L' C^ R'}", "lower_th_c"),[])

    def test_upper_thc_regex(self):
        self.assertEqual(regex_match("<L E^* R> ", "upper_th_c"),['E'])
        self.assertEqual(regex_match("<A B^ C>{1 2^* 3}<D E^* F>", "upper_th_c"),['E'])
        self.assertEqual(regex_match("{L' A^* R'}<1 2^* 3><M' C^* R'>", "upper_th_c"),['2','C'])
        self.assertEqual(regex_match("<1 2^ 3>{L' C^ R'}", "upper_th_c"),[])

    def test_empty(self):
        assert regex_match("", "upper_th_c") == []

if __name__ == '__main__':
    unittest.main()
