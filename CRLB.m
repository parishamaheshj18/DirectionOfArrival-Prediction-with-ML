clear all;
clc;
Data = csvread('CovTargetTrainingSNR2Range0.9NumSig3.csv');
Elements = 4;
NumSignals = 3;
d_lambda = 0.5;
mat = zeros(Elements-1,1);
SNR =2;
CRLB_value = zeros(3,3);
Mean_CRLB=zeros(5,1);
for SNR=0:2:10
%for k=1:length(Data)
for i=1: (Elements)
mat(i,1) = j.*2.*pi.*d_lambda.*(i-1);   
end

mat_1 = mat .*sin(Data(1,:));
exp_mat = exp(mat_1);
A = exp_mat;
sigma_sqr = (1/db2pow(SNR));

R = A*ctranspose(A) + sigma_sqr.*eye(Elements,Elements);

for i=1:Elements
 diff(i,:) = 1j*2*pi*d_lambda * cos(Data(1,:))*(i-1);
end

A_dot = A.*diff;

H = ctranspose(A_dot)  *  ( eye(Elements)- ( A *  inv(ctranspose(A)*A)  * ctranspose(A)) )  *A_dot;

CRLB_value = (sigma_sqr/(2*NumSignals)) * real( inv(H .* transpose(ctranspose(A) * inv(R) * A)) );
Mean_CRLB = sum(diag(CRLB_value));
plot(SNR,Mean_CRLB,'*r');
xlabel("SNR(db)");
ylabel("Mean CRLB");
hold on

end
