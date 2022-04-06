# test_dummy.py
'''
RSA 6/4/22
---------------------------
This demonstrates a test file that is run autimatically in integration. 
Any function thatstarts test_ in a file that startes test_ will be run automatically following a pull request.
---------------------------
'''

def capitalize_string(s):
    if not isinstance(s, str):
        raise TypeError('Please provide a string')
    return s.capitalize()

def test_capitalize_string():
    assert capitalize_string('test') == 'Test'
