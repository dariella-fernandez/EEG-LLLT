%%Check if tls and pbo baselines are statistically different

mean_tls_base = mean(pxx_tls_base, 3);
mean_pbo_base = mean (pxx_pbo_base, 3);
[h p] = ttest2(mean_tls_base',mean_pbo_base','alpha',0.05);
sorted_p = sort(p); 

channel = 34;
figure; 
plot(f1,10*log10(mean(pxx_pbo_base120(channel,:,:),3)), 'color','k')
hold on;
plot(f1,10*log10(mean(pxx_pbo_first(channel,:,:),3)), 'color','r')
hold on;
plot(f1,10*log10(mean(pxx_pbo_second(channel,:,:),3)), 'color','b')
hold on;
plot(f1,10*log10(mean(pxx_pbo_rec(channel,:,:),3)), 'color','m')
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
xlim([1 70])
title("Avg. PSD for PBO: all time periods")
legend('baseline','min 1-4', 'min 4-8', 'recovery')