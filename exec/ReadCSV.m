fileToRead1 = '/home/eds/research/somar/exec/TE_slice.csv';

% Import the file
newData1 = importdata(fileToRead1);

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
  assignin('base', vars{i}, newData1.(vars{i}));
end

% Create a variable for each column of data
for i = 1:length(colheaders)
  assignin('base', colheaders{i}, data(:,i));
end

% Free memory + clean up workspace
clear vars textdata newData1 fileToRead1 i data colheaders;

% Compute total energy
intTE = intKE + intAPEt;

% Normalize energies
normIntKE   = (  intKE(:) -   intKE(1)) / abs(  intKE(1));
normIntPE   = (  intPE(:) -   intPE(1)) / abs(  intPE(1));
normIntAPE  = ( intAPE(:) -  intAPE(1)) / abs( intAPE(1));
normIntAPEt = (intAPEt(:) - intAPEt(1)) / abs(intAPEt(1));
normIntTE   = (  intTE(:) -   intTE(1)) / abs(  intTE(1));

% Plot normalized energies
plot(time, normIntTE, '-o');
xlabel('time (s)');
ylabel('(E(t) - E(0)) / |E(0)|');