%function to check separability
cd ~/cvx;
cvx_setup;

%--->>>Parameter to check separation
eps = 1e-5;
%--->>>End parameter

%read data
X = dlmread(strcat(fileloc, 'Xs.txt'));
y = dlmread(strcat(fileloc,'Ys.txt'));
y1 = 2 * y -1;
p = size(X,2); %number of parameters

cvx_begin quiet
	variable betaHat(p)
	maximize( sum(y1 .* (X * betaHat)))
	subject to
		y1.*(X * betaHat) >= 0;
		-1 <= betaHat <= 1;
cvx_end

norm(betaHat)

if norm(betaHat)<eps
	sep = 0;	
else
	sep = 1;
end

filename = strcat(fileloc,'sep.txt');
fileID = fopen(filename,'w');
fprintf(fileID,'%d',sep);
fclose(fileID);

