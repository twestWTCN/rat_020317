ls = [-10:1:10];
for i = 1:length(ls) 

    rejrate(:,i) = sum((A>ls(i)),1)./size(A,1)';
    
end
figure
plot(ls,rejrate')
xlabel('Evidence Threshold 2ln(K)');
ylabel('Proportion of Exponents Valid');
legend(R.bandnames)