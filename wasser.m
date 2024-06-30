function [W2,dW2,dW2dt] = wasser(source_cdf,target_cdf,source_x,target_x)
% This function is based on the optimal-transport algorithm of Malcolm Sambridge 
% (Australian National University)
% Paper: 
% M. Sambridge, A. Jackson, and A. P. Valentine (2022), Geophysical inversion and optimal
% transport, Geophysical Journal International, 231(1), 172–198.
% doi: https://doi.org/10.1093/gji/ggac151
% Source code: 
% https://github.com/msambridge/waveform-ot 
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

a = [source_cdf(1:end-1);target_cdf];
[tk,tkarg] = sort(a);

indf = zeros(size(tk));
for in=1:length(tk)
    [~,tempel] = min(abs(source_cdf-tk(in)));
    if (abs(tk(in)-source_cdf(tempel))<10^(-12) && tk(in)>source_cdf(tempel))
        indf(in)=tempel-1;
    elseif tk(in)>source_cdf(tempel)
        indf(in)=tempel+1;
    else
        indf(in)=tempel;
    end
    if indf(in)<1 % This if-statement was added by jhunziker
        indf(in)=1;
    end
end

indg = zeros(size(tk));
for in=1:length(tk)
    [~,tempel] = min(abs(target_cdf-tk(in)));
    if (abs(tk(in)-target_cdf(tempel))<10^(-12) && tk(in)>target_cdf(tempel))
        indg(in)=tempel-1;
    elseif tk(in)>target_cdf(tempel)
        indg(in)=tempel+1;
    else
        indg(in)=tempel;
    end
    if indg(in)<1 % This if-statement was added by jhunziker
        indg(in)=1;
    end
end

dtk = [tk(1);tk(2:end) - tk(1:end-1)];

xft = source_x(indf);
xgt = target_x(indg);
dxft = abs(xft-xgt);

B = triu(ones(length(source_cdf),length(target_cdf)));
C = B-source_cdf.';
D = [C(:,1:end-1),zeros(length(source_cdf),length(target_cdf))];
Difftk = D(:,tkarg);
Diffdtk = [Difftk(:,1),Difftk(:,2:end)-Difftk(:,1:end-1)];

% Calculate W2^2 and its derivative with respect to unnormalised source amplitudes
dsqxft = dxft.*dxft;
W2 = sum(dsqxft.*dtk);
dsqxftdg = 2.0*(xft-xgt); % derivative of distances with respect to translation of input PDF position.
dW2 = sum(Diffdtk.*repmat(dsqxft.',[length(source_cdf),1]),2);
dW2dt = sum(dsqxftdg.*dtk);
