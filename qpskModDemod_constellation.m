%XXXX QPSK Modulation and Demodulation without consideration of noise XXXXX
clear all; close all;
data = round(rand(1,1000)); 
data2 = round(rand(1,1000));
% stem(data, 'linewidth',1.5), grid on; %stem:discrete plotting
% title('Transmiting signal'); %axis([ 0 20 0 1.5]);

data_NZR=2*data-1; 
spdata=reshape(data_NZR,2,length(data)/2); 
data_NZR2=2*data2-1; 
spdata2=reshape(data_NZR2,2,length(data2)/2);
br=10^6; %let us transmission bit rate 1000000
f=br; %minimum carrier frequency
T=1/br; %bit duration
t=T/99:T/99:T; %time vector for one bit information

% XXXXXXXXXXXXXXXXXXXXXXX QPSK modulation XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
y=[]; y_in=[]; y_qd=[]; yy=[]; y_in2=[]; y_qd2=[];
const = [];
for ii=1:length(data)/2
    % 每個data分別*cos(2pi)的vector,cos(2pi*[1/99,2,99...99/99])
    % ft=[1/99,2/99...99/99]
    y1=spdata(1,ii)*cos(2*pi*f*t); %inphase component
    y2=spdata(2,ii)*sin(2*pi*f*t) ;%quadrature component
    y=[y y1+y2]; 
    if spdata(1,ii)==-1 && spdata(2,ii)==-1 %point 00
        C= exp(((pi)/4)*1i);
    elseif spdata(1,ii)==-1 && spdata(2,ii)==1 %point 01
        C = exp(((3*pi)/4)*1i);
    elseif spdata(1,ii)==1 && spdata(2,ii)==-1 %point 10
        C = exp(((5*pi)/4)*1i);
    elseif spdata(1,ii)==1 && spdata(2,ii)==1 %point 11
        C = exp(((7*pi)/4)*1i);
    end
    const = [const C];
end
Tx_sig=y;
for i=1:length(data2)/2
    y12=spdata2(1,i)*cos(2*pi*f*t); 
    y22=spdata2(2,i)*sin(2*pi*f*t) ;
    y_in2=[y_in2 y12]; 
    y_qd2=[y_qd2 y22]; 
    yy=[yy y12+y22];
end
Tx_sig2=yy;

tt=T/99:T/99:(T*length(data))/2;
snr=0.01;
txsig = Tx_sig;
Tx_sig = awgn(Tx_sig,snr,'measured','linear');
Tx_sig2 = awgn(Tx_sig2,snr,'measured','linear');
% noise = Tx_sig - txsig;
% pow_noise = mean(noise.^2)

% figure(1)
% subplot(3,1,1);
% plot(tt,y_in,'linewidth',1.5), grid on;
% title(' Inphase component in QPSK modulation ');
% xlabel('time(sec)'); ylabel(' amplitude(volt)');
% subplot(3,1,2);
% plot(tt,y_qd,'linewidth',1.5), grid on;
% title(' Quadrature component in QPSK modulation ');
% xlabel('time(sec)'); ylabel(' amplitude(volt)');
% subplot(3,1,3);
% plot(tt,Tx_sig,'r','linewidth',1.5), grid on;
% title('QPSK modulated signal (sum of inphase and Quadrature phase signal)');
% xlabel('time(sec)'); ylabel(' amplitude(volt)');
% figure(2)
% subplot(3,1,1);
% plot(tt,y_in2,'linewidth',3), grid on;
% title(' 2 Inphase component in QPSK modulation ');
% xlabel('time(sec)'); ylabel('amplitude(volt)');
% subplot(3,1,2);
% plot(tt,y_qd2,'linewidth',3), grid on;
% title(' 2 Quadrature component in QPSK modulation ');
% xlabel('time(sec)'); ylabel('amplitude(volt)');
% subplot(3,1,3);
% plot(tt,Tx_sig2,'r','linewidth',3), grid on;
% title('2 QPSK modulated signal (sum of inphase and Quadrature phase signal)');
% xlabel('time(sec)'); ylabel('amplitude(volt)');


% XXXXXXXXXXXXXXXXXXXXXXXXXXXX ANTENNA ARRAY XXXXXXXXXXXXXXXXXXXXXXXXXX
dis=0.5; num=10; arr=0:num-1;  %0~9
deg = [-50,20]; rad=deg/180*pi;
a1 = exp(-2j*pi * dis * arr * sin(rad(1))); %一個角度對每個天線
a2 = exp(-2j*pi * dis * arr * sin(rad(2)));
a=[a1;a2];
r = a' * [Tx_sig;Tx_sig2];
% XXXXXXXXXXXXXXXX
thetas=-1/2*pi:pi/300:1/2*pi; results = [];
for i=thetas
    w = exp(-2j*pi*dis*arr*sin(i));
    r_weighted = w * r;
    results=[results mean(abs(r_weighted))];
end 
% plot the distribution of results
% plot(thetas*180/pi,results);xlabel('Direction (degree)');ylabel('Magnitude');
%find max
[~,arg]=max(results);
g=thetas(arg)*180/pi;
%find 2nd max (suppose out of ±10deg)
for i=2:1:length(results)
    [~,argg]=maxk(results,i);
    gg=thetas(argg(i))*180/pi;
    if ((g-10)>gg | gg>(g+10))==1
        break
    end
end
dir=[g gg];
% XXXXXXXXXXXXXXXX
a1h = ctranspose(a1);
m1 = a1h.' * a1.';
a2h = ctranspose(a2);
m2 = a2h.' * a2.';
% array receiving sig
r = a1.'*Tx_sig;
r = r + a2.'*Tx_sig2; %如果第二個direction是傳進第二個source 這裡要改
%reconstruct tx_sig≈s1
s1 = a1h.' * r /m1; s1=real(s1);
s2 = a2h.' * r /m2; s2=real(s2);

% % XXXXXXXXXXXXXXXXXXXXXXXXXXXX QPSK demodulation XXXXXXXXXXXXXXXXXXXXXXXXXX
Rx_data=[];
Rx_sig=s1;
arrr=[]; const2=[];
for ii=1:1:length(data)/2
    %%XXXXXX inphase coherent dector XXXXXXX
    % Rx_sig shape(1,99*30) % for每個i 照順序取出Rx_sig裡面的99個time vector
    % (.^對應元素相乘)乘出來是(1,99)的vector
    Z_in=Rx_sig((ii-1)*length(t)+1:ii*length(t)).*cos(2*pi*f*t);
%     Z_in_intg=((T/99)*trapz(Z_in))*(2/T); %integration using trapizodial rule
    Z_in_intg=(trapz(Z_in)/198);
    arrr=[arrr Z_in_intg];
    if(Z_in_intg>0) 
        Rx_in_data=1;
    else
        Rx_in_data=0;
    end    
    %%XXXXXX Quadrature coherent dector XXXXXX
    Z_qd=Rx_sig((ii-1)*length(t)+1:ii*length(t)).*sin(2*pi*f*t);
    Z_qd_intg=(trapz(Z_qd)/198);
    const2 = [const2 Z_qd_intg];
    if (Z_qd_intg>0)
        Rx_qd_data=1;
    else
        Rx_qd_data=0; 
    end
    Rx_data=[Rx_data  Rx_in_data  Rx_qd_data];
    if Rx_in_data==0 && Rx_qd_data==0 %point 00
        C= exp(1i*((pi)/4));
    elseif Rx_in_data==0 && Rx_qd_data==1 %point 01
        C = exp(1i*((3*pi)/4));
    elseif Rx_in_data==1 && Rx_qd_data==0 %point 10
        C = exp(1i*((5*pi)/4));
    elseif Rx_in_data==1 && Rx_qd_data==1 %point 11
        C = exp(1i*((7*pi)/4));
    end
%     const2 = [const2 C];
end
const = awgn(const,24,'measured');
scatterplot(const); grid minor;
% const2 = awgn(const2,20,'measured');
% scatterplot(const2); grid minor;

Rx_data2=[];
Rx_sig2=s2;
for i=1:1:length(data2)/2
    Z_in2=Rx_sig2((i-1)*length(t)+1:i*length(t)).*cos(2*pi*f*t); 
    Z_in_intg2=(trapz(t,Z_in2))*(2/T);
    if(Z_in_intg2>0)
        Rx_in_data2=1;
    else
       Rx_in_data2=0;
    end
    Z_qd2=Rx_sig2((i-1)*length(t)+1:i*length(t)).*sin(2*pi*f*t);
    Z_qd_intg2=(trapz(t,Z_qd2))*(2/T);
    if (Z_qd_intg2>0)
        Rx_qd_data2=1;
    else
        Rx_qd_data2=0; 
    end
    Rx_data2=[Rx_data2  Rx_in_data2  Rx_qd_data2];
end
% figure(1)
% subplot(211)
% stem(data, 'linewidth',1.5), grid on; %stem:discrete plotting
% title('data (Transmiting signal)'); %axis([ 0 20 0 1.5]);
% subplot(212)
% stem(Rx_data,'linewidth',1.5) 
% title('Rx_data (Receiveing signal)'); %axis([ 0 20 0 1.5]), grid on;
% % figure(2)
% % subplot(211)
% % stem(data2, 'linewidth',3), grid on;
% % title('2 Transmiting signal'); axis([ 0 20 0 1.5]);
% % subplot(212)
% % stem(Rx_data2,'linewidth',3) 
% % title('2 Receiveing signal'); axis([ 0 20 0 1.5]), grid on;
% 
% figure(1)
% subplot(211)
% plot(tt,Tx_sig,'linewidth',1), grid on;
% title('tx'); xlabel('t'); ylabel(' amplitude'); %set(gca, 'ylim', [-3,3]);
% subplot(212)
% plot(tt,s1,'linewidth',1), grid on;
% title('s1'); xlabel('t'); ylabel(' amplitude');
% figure(2)
% subplot(211)
% plot(tt,Tx_sig2,'linewidth',2), grid on;

% title('tx2');
% subplot(212)
% plot(tt,s2,'linewidth',2), grid on; 
% title('s2');  

[xx,yy] = biterr(data,Rx_data);
