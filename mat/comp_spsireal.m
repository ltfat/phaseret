function [chat,endphase] = comp_spsireal(s,a,M,startphase)
M2 = floor(M/2) + 1;
N = size(s,2);

if nargin<4 || isempty(startphase)
    m_phase = zeros(size(s,1),1); 
else
    m_phase = startphase;
end

chat = zeros(size(s));

for n=1:N
    sabs=abs(s(:,n));    %magnitude of the fft
    centno = 0;
    %   check through the magnitude spectrum to find peak
    for m=2:1:M2-1
        % If the current bin is bigger than both his neighbors
        if(sabs(m)>sabs(m-1) && sabs(m)>sabs(m+1))
            alpha=log(sabs(m-1)+eps);
            beta=log(sabs(m)+eps);
            gamma=log(sabs(m+1)+eps);
            % This is 2nd difference
            denom=alpha-2*beta+gamma;

            % initialize p
            if(denom~=0)
                % Peak offset from m
                p=0.5*(alpha-gamma)/denom;
            else
                % denom is zero only at inflection point (2nd deriv is zero)
                p=0;
            end
            
            % Peak is at m+p in Matlab indexing
            instf = (m-1+p);
            centno = centno + 1;
            peakPhase = m_phase(m) + 2*pi*a*instf/M; %phase accumulator for this peak bin
            % Assign the peak phase to current bin
            m_phase(m) = peakPhase;

          
            % Go down from the peak until valey
            bin=m-1;
            while(bin>1 && sabs(bin)<sabs(bin+1))
                m_phase(bin)=peakPhase;
                bin=bin-1;
            end
              
            % Go up from the peak until valey
            bin=m+1;
            while(bin<(M2) && sabs(bin)<sabs(bin-1))
                m_phase(bin)=peakPhase;
                bin=bin+1;
            end
            
        end
    end
    
chat(:,n) = sabs.*exp(1i*m_phase);    
end

endphase = m_phase;
