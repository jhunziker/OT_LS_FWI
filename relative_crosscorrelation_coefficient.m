function Tfrac = relative_crosscorrelation_coefficient(wobs,wpred,tobs,srcplot,recplot)
% Calculate the relative crosscorrelation criterion. 
%
% This file is part of the OT-LS-FWI package. OT-LS-FWI is free software: 
% you can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. 
% OT-LS-FWI is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details. You should have received a copy of the GNU General Public 
% License along with OT-LS-FWI. If not, see <https://www.gnu.org/licenses/>. 

fs = 12; % Fontsize
lw = 2; % Linewidth

Tfrac = zeros(size(wobs,2),size(wobs,3));
for srcnum=1:size(wobs,3)
    for recnum=1:size(wobs,2)

        % DETERMINE PERIOD OF PREDICTED TRACE
        % (If noise is present, the local maxima are only correctly determined, 
        %  if there is no noise, hence in the predicted data and not in the observed data.)
        % Select trace
        wave = squeeze(wpred(:,recnum,srcnum));
        % Find the local maxima
        locmaxvec = islocalmax(wave);
        % Find the largest local maxima
        wavelocmax = wave(locmaxvec);
        tobslocmax = tobs(locmaxvec);
        [~,maxel] = max(wavelocmax);
        % Check if the left or right local maxima is larger
        if maxel==1
            maxright = wavelocmax(maxel+1);
            t1 = tobslocmax(maxel);
            a1 = wavelocmax(maxel); % only for plotting
            t2 = tobslocmax(maxel+1);
            a2 = wavelocmax(maxel+1); % only for plotting
	elseif maxel==length(wavelocmax)
	    maxleft = wavelocmax(maxel-1);
            t1 = tobslocmax(maxel-1);
            a1 = wavelocmax(maxel-1); % only for plotting
            t2 = tobslocmax(maxel);
            a2 = wavelocmax(maxel); % only for plotting
        else
	    maxleft = wavelocmax(maxel-1);
            maxright = wavelocmax(maxel+1);
            if wavelocmax(maxel-1)>wavelocmax(maxel+1)
                t1 = tobslocmax(maxel-1);
                a1 = wavelocmax(maxel-1); % only for plotting
                t2 = tobslocmax(maxel);
                a2 = wavelocmax(maxel); % only for plotting
            else
                t1 = tobslocmax(maxel);
                a1 = wavelocmax(maxel); % only for plotting
                t2 = tobslocmax(maxel+1);
                a2 = wavelocmax(maxel+1); % only for plotting
            end
        end
        % Determine the period in seconds
        period = t2-t1;

        if (srcnum==srcplot && recnum==recplot)
            figure;
            subplot(311)
            plot(tobs*1e9,wobs(:,recnum,srcnum),'k','Linewidth',lw)
            hold on
            plot(tobs*1e9,wpred(:,recnum,srcnum),'r--','Linewidth',lw)
            % plot(tobs(locmaxvec)*1e9,wobs(locmaxvec,recnum,srcnum),'ko')
            plot(t1*1e9,a1,'ko');
            plot(t2*1e9,a2,'ko');
            hold off
            xlabel('Time (ns)')
            title(['Source ',num2str(srcnum),'; Receiver ',(num2str(recnum))])
            set(gca,'Fontsize',fs)
        end

        % Do the crosscorrelation
        corrvec = xcorr(squeeze(wobs(:,recnum,srcnum)),squeeze(wpred(:,recnum,srcnum)));
        % Determine how many elements does the maximum deviate from the center
        [~,maxcorrel] = max(corrvec);
        shiftel = maxcorrel-(length(corrvec)-1)/2;
        if (srcnum==srcplot && recnum==recplot)
            % Plot crosscorrelation
            subplot(312)
            plot(corrvec,'b','Linewidth',lw)
            xlim([1,length(corrvec)])
            set(gca,'Fontsize',fs)

            % Shift one trace accordingly
            subplot(313)
            plot(tobs*1e9,wobs(:,recnum,srcnum),'k','Linewidth',lw)
            hold on
            plot(tobs*1e9,circshift(squeeze(wpred(:,recnum,srcnum)),shiftel),'r--','Linewidth',lw)
            hold off
            xlabel('Time (ns)')
            set(gca,'Fontsize',fs)
        end

        Tfrac(recnum,srcnum) = shiftel*(tobs(2)-tobs(1))/period;
    end
end
