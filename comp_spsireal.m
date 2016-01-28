function [chat,newphase] = comp_spsireal(c,a,M,mask,startphase)
M2 = floor(M/2) + 1;
N = size(c,2);

if nargin <5 || isempty(startphase)
    m_phase = zeros(size(c,1),1); 
else
    m_phase = startphase;
end


newphase = zeros(size(c));

for n=1:N
    sabs=abs(c(:,n));    %magnitude of the fft
    centno = 0;
%   check through the magnitude spectrum to find peak
    for m=2:1:M2-1
        if(sabs(m)>sabs(m-1) && sabs(m)>sabs(m+1))

        %  **************** ORIGINAL COMMENT *****************************
        %  Phase estimate using quadratic interpolation as per Julios O. Smith
        %  http://www.dsprelated.com/dspbooks/sasp/Quadratic_Interpolation_Spectral_Peaks.html
        %  Use quadratic interpolation to estimate the real peak position
        %  ***************************************************************
            alpha=log(sabs(m-1)+eps);
            beta=log(sabs(m)+eps);
            gamma=log(sabs(m+1)+eps);
            denom=alpha-2*beta+gamma;

            % initialize p
            if(denom~=0)
                p=0.5*(alpha-gamma)/denom;
            else
                p=0;
            end
            
            instf = (m-1+p);
            centno = centno + 1;
            peakPhase = m_phase(m) + 2*pi*a*instf/M; %phase accumulator for this peak bin
            m_phase(m)=peakPhase;

            %  **************** ORIGINAL COMMENT *****************************
            % Apply_input simple phase locking around the peaks.
            % The phase relationships of the bins around the peak were determined by
            % some simple experiments inspired by: 
            % Laroche/Dolson "About This Phasiness Business" (1997), which mentions a paper
            % M.S. Puckette "Phase-locked vocoder" (1995).  
            % http://msp.ucsd.edu/Publications/mohonk95.pdf

            % According to Laroche/Dolson:
            % "Puckette in [5] recognized that for a constant-frequency_input constant-amplitude sinusoid 
            % the sy_inputnthesis phases around the maximum of the Fourier transform should exhibit +/- pi 
            % alternations and proposed a very_input simple way_input to constrain them to do so".
            %  ***************************************************************
            
            %  **************** NEW COMMENT *******************************
            % This is all not relevant when computing dgt properly!
            % Therefore there is now +0 everywhere

        % If actual peak is to the right
            if (p>0)

            % First bin to right has pi shift
                bin=m+1;
                m_phase(bin)=peakPhase+0;
                bin=m-1;

                while(bin>1 && sabs(bin)<sabs(bin+1))
                    % Bins to left have shift of pi
                    m_phase(bin)=peakPhase+0;
                    bin=bin-1;
                end
                bin=m+2;

                while(bin<(M2) && sabs(bin)<sabs(bin-1))
                    % Other bins (apart from the immediate bin to the right of the peak)
                    % have zero shift
                    m_phase(bin)=peakPhase+0;
                    bin=bin+1;
                end
            end

            % Peak is to the left
             if(p<0)
            % First bin to left has pi shift
                bin=m-1;
                m_phase(bin)=peakPhase+0;

                bin=m+1;
                while(bin<(M2) && sabs(bin)<sabs(bin-1))
                    % Bins to right have shift of pi
                    m_phase(bin)=peakPhase+0;
                    bin=bin+1;
                end
                 bin=m-2;
                 while(bin>1 && sabs(bin)<sabs(bin+1))
                    % Other bins to left have zero shift
                    m_phase(bin)=peakPhase+0;
                    bin=bin-1;
                 end
             end
        end
    end
    
   if ~isempty(mask) 
       % Use the reliable phase
       m_phase(mask(:,n)) = angle(c(mask(:,n),n));
   end
   newphase(:,n) = m_phase; 
end

chat = abs(c).*exp(1i*newphase);