channel = 34;
subject = 2;
figure;
plot(f1,10*log10(pxx_tls_base(channel,:,1)), 'color','r')
hold on;
plot(f1,10*log10(pxx_tls_base(channel,:,2)), 'color','g')
hold on;
plot(f1,10*log10(pxx_tls_base(channel,:,3)), 'color','b')
hold on;
plot(f1,10*log10(pxx_tls_base(channel,:,4)), 'color','m')
hold on;
plot(f1,10*log10(pxx_tls_base(channel,:,5)), 'color','c')
hold on;
plot(f1,10*log10(pxx_tls_base(channel,:,6)), 'color','y')
hold on
plot(f1,10*log10(pxx_tls_base(channel,:,7)), 'color','k')
hold on;
plot(f1,10*log10(pxx_tls_base(channel,:,8)), 'color','r')
hold on;
plot(f1,10*log10(pxx_tls_base(channel,:,9)), 'color','g')
hold on;
plot(f1,10*log10(pxx_tls_base(channel,:,10)), 'color','b')
hold on;
plot(f1,10*log10(pxx_tls_base(channel,:,11)), 'color','m')
hold on;
plot(f1,10*log10(pxx_tls_base(channel,:,12)), 'color','c')
hold on;
plot(f1,10*log10(pxx_tls_base(channel,:,13)), 'color','y')
hold on;
plot(f1,10*log10(pxx_tls_base(channel,:,14)), 'color','k')
hold on
plot(f1,10*log10(pxx_tls_base(channel,:,15)), 'color','k')
hold on;
%plot(f1,10*log10(pxx_tls_base(channel,:,16)), 'color','g')
hold on;
plot(f1,10*log10(pxx_tls_base(channel,:,17)), 'color','g')
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
xlim([1 70])
title("Baseline PSD for all TLS subjects using last 60 secs")
% legend('rec','baseline')