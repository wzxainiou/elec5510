clc
clear all
close all
tic

%Elapsed time is  seconds.

N_iter=100; % number of iterations for each SNR
n_ldpc=64800; % codeword length of LDPC
R=[2/5,9/10,2/3,8/9]; % LDPC code rate

figure(1)

for codecase=1:4
    
    Rc=R(codecase);
    k_ldpc=n_ldpc*Rc;

    switch codecase
    case 1
        Eb_N0=[(-4:1:-1),(-0.9:0.1:0),(0.1:0.05:0.4),(0.41:0.01:0.5)];
        M=4;
    case 2
        Eb_N0=[(-4:1:3),(3.1:0.1:3.6),(3.61:0.01:3.7),(3.705:0.005:3.8)];
        M=4;
    case 3
        Eb_N0=[(-4:1:3),(3.1:0.1:3.5),(3.51:0.01:4)];
        M=8;
    case 4
        Eb_N0=[(-4:1:5),(5.2:0.2:6.2),(6.21:0.01:6.34),(6.345:0.005:6.400)];
        M=8;
    end
    
    fprintf("Case %d: LDPC code rate =%s %dPSK",codecase,rats(Rc),M);
    
    % LDPC Encoding and Decoding Objects
    hEnc = comm.LDPCEncoder('ParityCheckMatrix',dvbs2ldpc(Rc));
    hDec = comm.LDPCDecoder('ParityCheckMatrix',dvbs2ldpc(Rc));
    
    % Modulation Object
    hMod = comm.PSKModulator(M, 'BitInput',true);
    
    n=1;
    BER=zeros(size(Eb_N0));
    
    for eb_no=Eb_N0
        
        [Rc, eb_no, M]
        
        gamma=eb_no+10*log10(log2(M)*Rc);
        
        % AWGN
        hChan = comm.AWGNChannel(...
            'NoiseMethod','Signal to noise ratio (SNR)','SNR',gamma);
        
        % Demodulation
        hDemod = comm.PSKDemodulator(M, 'BitOutput',true,...
            'DecisionMethod','Approximate log-likelihood ratio', ...
            'Variance', 1/10^(hChan.SNR/10));
        
        % For BER calculations
        hError = comm.ErrorRate;
        
        EE=zeros(1,N_iter);
        
        for counter = 1:N_iter
            reset(hEnc);
            reset(hDec);
            reset(hError);
            
            data = logical(randi([0,1],k_ldpc, 1)); % Random bit generation 
            encodedData = step(hEnc,data); % LDPC encoding
            modSignal = step(hMod, encodedData); % Modulation
            receivedSignal = step(hChan, modSignal); % Transmission over the channel
            demodSignal = step(hDemod, receivedSignal); % Demodulation
            RBits = step(hDec, demodSignal); % LDPC decoding
            errorStats = step(hError, data, RBits); % BER calculations
            EE(1,counter)=errorStats(1);
        end
        BER(1,n)=mean(EE(1,:)); % Averaging over iterations
        n=n+1;
    end
    
    semilogy(Eb_N0,BER)
    hold on

end

%Case 5&6:Uncoded
for codecase=5:6
    
    switch codecase
        case 5
            Eb_N0=(-4:1:11);
            M=4;
        case  6
            Eb_N0=(-4:1:11);
            M=8;
    end
    
    fprintf("Case %d:Uncoded, %dPSK",codecase,M);
    
    % Modulation Object
    hMod = comm.PSKModulator(M, 'BitInput',true);
    
    n=1;
    BER=zeros(size(Eb_N0));
    
    for eb_no=Eb_N0
        
        [1, gamma, M]
        
        gamma=eb_no+10*log10(log2(M));
        
        hChan = comm.AWGNChannel(...
            'NoiseMethod','Signal to noise ratio (SNR)','SNR',gamma);
        hDemod = comm.PSKDemodulator(M, 'BitOutput',true,...
            'DecisionMethod','Approximate log-likelihood ratio', ...
            'Variance', 1/10^(hChan.SNR/10));
        hError = comm.ErrorRate; % For BER calculations
        EE=zeros(1,N_iter);
        
        for counter = 1:N_iter
            
            reset(hError);
            
            data = logical(randi([0,1],64800, 1)); % Random bit generation
            encodedData = data;
            
            modSignal = step(hMod, encodedData); % Modulation
            receivedSignal = step(hChan, modSignal); % Transmission over the channel
            demodSignal = step(hDemod, receivedSignal); % Demodulation
            
            hard_symbols = (demodSignal > 0) + (demodSignal <= 0).*(-1); % Hard Decoding without LDPC
            demodSignal = (hard_symbols == -1);
            receivedBits = demodSignal;
            
            errorStats = step(hError, data, receivedBits); % BER calculations
            EE(1,counter)=errorStats(1);
        end
        
        BER(1,n)=mean(EE(1,:)); % Averaging over iterations
        n=n+1;
        
    end
    
    semilogy(Eb_N0,BER)
    hold on
end 

grid on;
title('LDPC simulation'); xlabel('Eb/N_0 [dB]'); ylabel('BER');
legend('R=2/5,QPSK','R=9/10,QPSK','R=2/3,8PSK','R=8/9,8PSK','Uncoded, QPSK','Uncoded, 8PSK');

toc
