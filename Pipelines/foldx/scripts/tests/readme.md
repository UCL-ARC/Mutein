### Continuous Integration and Testing

- CI tests are in folders starting tests_ and are functions starting tests_
- Manually run tests are in folders manually_ these should be run before a checking
- Tests that capture a particular case to reproduce are otherwise named 

The entire foldx application cannot be run, or other external libraries, in CO on github (????) but we can assume they will remain stable if the function call parameters remain the same, so we test the inputs to these functions do not change in CI.