Instructions to Run the Matlab code:

1) Generating a new instance:
To create a new instance run the following command

<InstanceName> = getStructure(<n>,<k>);

for generating dense or sparse matrix you only have to comment one of the lines 21 and 24 in genStructure.m file in sparse case you have to specify the value of desired density.

2) Starting the solver PDIPM:
To start optimization by PDIPM run the following command

[primal, x, status] = PDIPM(<InstanceName>,<MaxIter>)


