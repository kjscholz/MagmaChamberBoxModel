function Trt = heat_conduction_chamber_profile(maxn,a,c,r, kappa,Tb)


global storeTime storeTemp storeSumk storeSumk_2 storeSumk_old storeSumk_2_old lengthTime maxTime
global switch_Tprofile


% temperature history
time = storeTime;
Tr0   = storeTemp;

time_index= length(time);

% pick a time
current_time    = time(end);

if time_index == 2
    % first sum over n
    
    sumn = 0;
    for n=1:maxn
        % sum over k within first sum over n
        sumk=0;
        for k=1:time_index-1
            past_time       = time(k);
            past_time_delta = time(k+1);
            
            mk       = (Tr0(k+1) - Tr0(k))/(time(k+1)-time(k));
            pk       = Tr0(k) - mk*time(k);
            
            termk = mk*c^4/(kappa^2*n^4*pi^4) ...
                * (exp(-kappa*(n*pi/c)^2*(current_time-past_time))*(1-kappa*(n*pi/c)^2*past_time) ...
                + exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta))*(kappa*(n*pi/c)^2*past_time_delta-1)) ...
                + pk*c^2/(kappa*n^2*pi^2) ...
                * (exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) ...
                - exp(-kappa*(n*pi/c)^2*(current_time-past_time)));
            sumk = sumk + termk;
        end
        termn = n.*sin(n.*pi.*(r-a)./c).*sumk;
        sumn  = sumn + termn;
    end
    
    
    term1 = -4.*pi.*kappa.*a./(2.*r.*c^2).*(-1).*sumn;
    
    
    % second sum over n
    sumn = 0;
    for n=1:maxn
        % sum over k within second sum over n
        sumk = 0;
        for k=1:time_index-1
            past_time       = time(k);
            past_time_delta = time(k+1);
            termk = exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) ...
                - exp(-kappa*(n*pi/c)^2*(current_time-past_time));
            sumk = sumk + termk;
        end
        termn = 1./n.*sin(n.*pi.*(r-a)./c).*cos(n*pi).*sumk;
        sumn  = sumn + termn;
    end
    
    
    term2 = -4.*pi.*kappa.*(a+c)./(2.*kappa.*pi^2.*r).*Tb.*(1).*sumn;
elseif time_index > 2 && time_index>lengthTime
    % first sum over n
    sumn = 0;
    for n=1:maxn
        % sum over k within first sum over n
        
        temp=storeSumk(n)*exp(-kappa*(n*pi/c)^2*(current_time-time(time_index-1)));
        
        past_time       = time(time_index-1);
        past_time_delta = time(time_index);
        
        mk       = (Tr0(time_index) - Tr0(time_index-1))/(time(time_index)-time(time_index-1));
        pk       = Tr0(time_index-1) - mk*time(time_index-1);
        
        termk = mk*c^4/(kappa^2*n^4*pi^4) ...
            * (exp(-kappa*(n*pi/c)^2*(current_time-past_time))*(1-kappa*(n*pi/c)^2*past_time) ...
            + exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta))*(kappa*(n*pi/c)^2*past_time_delta-1)) ...
            + pk*c^2/(kappa*n^2*pi^2) ...
            * (exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) ...
            - exp(-kappa*(n*pi/c)^2*(current_time-past_time)));
        sumk = temp + termk;

        termn = n.*sin(n.*pi.*(r-a)./c).*sumk;
        sumn  = sumn + termn;
    end
    
    
    term1 = -4.*pi.*kappa.*a./(2.*r.*c^2).*(-1).*sumn;

    % second sum over n
    sumn = 0;
    for n=1:maxn
        temp=storeSumk_2(n)*exp(-kappa*(n*pi/c)^2*(current_time-time(time_index-1)));
        past_time       = time(time_index-1);
        past_time_delta = time(time_index);
        termk = exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) ...
            - exp(-kappa*(n*pi/c)^2*(current_time-past_time));
        sumk = temp + termk;

        termn = 1./n.*sin(n.*pi.*(r-a)./c).*cos(n*pi).*sumk;
        sumn  = sumn + termn;
    end
    term2 = -4.*pi.*kappa.*(a+c)./(2.*kappa.*pi^2.*r).*Tb.*(1).*sumn;
        switch_Tprofile=0;
elseif time_index > 2 && time_index==lengthTime
    sumn = 0;

      for n=1:maxn
        % sum over k within first sum over n
        if switch_Tprofile==0
            temp=storeSumk_old(n)*exp(-kappa*(n*pi/c)^2*(current_time-time(end-1)));
            past_time       = time(time_index-1);
            past_time_delta = time(time_index);
            
            mk       = (Tr0(time_index) - Tr0(time_index-1))/(time(time_index)-time(time_index-1));
            pk       = Tr0(time_index-1) - mk*time(time_index-1);
            
            termk = mk*c^4/(kappa^2*n^4*pi^4) ...
                * (exp(-kappa*(n*pi/c)^2*(current_time-past_time))*(1-kappa*(n*pi/c)^2*past_time) ...
                + exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta))*(kappa*(n*pi/c)^2*past_time_delta-1)) ...
                + pk*c^2/(kappa*n^2*pi^2) ...
                * (exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) ...
                - exp(-kappa*(n*pi/c)^2*(current_time-past_time)));
            sumk = temp + termk;
        else
            sumk=storeSumk_old(n);
        end
        
        termn = n.*sin(n.*pi.*(r-a)./c).*sumk;
        sumn  = sumn + termn;
    end
    
    term1 = -4.*pi.*kappa.*a./(2.*r.*c^2).*(-1).*sumn;
    
    % second sum over n
    sumn = 0;
    for n=1:maxn
        if switch_Tprofile==0
            
            temp=storeSumk_2_old(n)*exp(-kappa*(n*pi/c)^2*(current_time-time(end-1)));
            past_time       = time(end-1);
            past_time_delta = time(end);
            termk = exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) ...
                - exp(-kappa*(n*pi/c)^2*(current_time-past_time));
            sumk = temp + termk;
            
        else
            sumk=storeSumk_2_old(n);
        end
        termn = 1./n.*sin(n.*pi.*(r-a)./c).*cos(n*pi).*sumk;
        sumn  = sumn + termn;
    end
    term2 = -4.*pi.*kappa.*(a+c)./(2.*kappa.*pi^2.*r).*Tb.*(1).*sumn;

elseif time_index > 2 && time_index<lengthTime
    switch_Tprofile=1;
    sumn = 0;
    %display([num2str(time_index)  ' ; ' num2str(time(end-2)) ' ; ' num2str(time(end-1)) ' ; ' num2str(time(end))]);
    for n=1:maxn        
        termn = n.*sin(n.*pi.*(r-a)./c).*storeSumk_old(n);
        sumn  = sumn + termn;
    end    
    
    term1 = -4.*pi.*kappa.*a./(2.*r.*c^2).*(-1).*sumn;

    % second sum over n
    sumn = 0;
    for n=1:maxn
        
        termn = 1./n.*sin(n.*pi.*(r-a)./c).*cos(n*pi).*storeSumk_2_old(n);
        sumn  = sumn + termn;
    end
    term2 = -4.*pi.*kappa.*(a+c)./(2.*kappa.*pi^2.*r).*Tb.*(1).*sumn;


elseif time_index <2
    term1 =zeros(size(r));
    term2 =zeros(size(r));
end

% third sum over n
Ta = Tr0(time_index);
sumn = 0;
for n=1:maxn
    
    termn =  sin(n.*pi.*(r-a)./c).*exp(-kappa.*(n.*pi./c)^2.*current_time) ...
        * ((a*(a+c)*Ta - a*(a+c)*Tb)/(n*pi)*(1-cos(n*pi)) ... %SIGN!!
        + ((a+c)*Tb -a*Ta)/(n*pi)*(-(a+c)*cos(n*pi) + a));
    
    sumn = sumn + termn;
end
term3 = 2./(r.*c).*sumn;

Trt = term1 + term2 + term3;

