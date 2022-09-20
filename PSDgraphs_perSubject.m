% % Individual Subject PSD Plots

channel = 34;
for sub=1:17
    figure;
    plot(f1,10*log10(pxx_tls_base120(channel,:,sub)), 'color','k')
    hold on;    
    plot(f1,10*log10(pxx_tls_first(channel,:,sub)), 'color','r')
    hold on;
    plot(f1,10*log10(pxx_tls_second(channel,:,sub)), 'color','b')
    hold on;
    plot(f1,10*log10(pxx_tls_rec(channel,:,sub)), 'color','m')
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
    xlim([1 70])
    title("PSD for TLS subject " + sub + " using last 60 secs")
    legend('baseline','min 1-4', 'min 4-8', 'recovery')
end

for sub=1:15
    figure;
    plot(f1,10*log10(pxx_pbo_base120(channel,:,sub)), 'color','k')
    hold on;    
    plot(f1,10*log10(pxx_pbo_first(channel,:,sub)), 'color','r')
    hold on;
    plot(f1,10*log10(pxx_pbo_second(channel,:,sub)), 'color','b')
    hold on;
    plot(f1,10*log10(pxx_pbo_rec(channel,:,sub)), 'color','m')
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
    xlim([1 70])
    title("PSD for PBO subject " + sub + " using last 60 secs")
    legend('baseline','min 1-4', 'min 4-8', 'recovery')
end