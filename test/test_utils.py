from autodeer.utils import gcd

def test_gcd():
    assert gcd([2.0, 4.0, 6.0, 8.0]) == 2
    assert gcd([1.5, 3.0]) == 1.5
    assert gcd([3, 6, 9, 12]) == 3
    assert gcd([5, 10, 15, 20]) == 5
    assert gcd([7, 14, 21, 28]) == 7
    assert gcd([8, 12, 16, 20]) == 4
    assert gcd([10, 20, 30, 40]) == 10
    assert gcd([15, 25, 35, 45]) == 5
    assert gcd([18, 24, 30, 36]) == 6
    assert gcd([21, 28, 35, 42]) == 7
    assert gcd([27, 36, 45, 54]) == 9