function [] = plot_onwarp(warpnoc, band, warptest,tg, s_sBn, plotallFlag, axes1)


% Parameter checking
if nargin < 2, errodlg('The number of parameters introduced are not expected.');end
if nargin < 3, warptest = []; end   
if nargin < 4 || isempty(tg), tg = 'x';end
if nargin < 5 || isempty(s_sBn), s_sBn = size(warpnoc,1);end
if nargin < 6,  plotallFlag = true;end
if nargin < 7, 
    h = figure;
    axes1 = axes; 
end;

s = size(warpnoc);
samplepoints = [1:s(1)]';
warpnocnew = warpnoc - repmat(samplepoints,1,s(2));
newband = zeros(s(1),2);

for i=1:s(1)
    newband(i,1) = band(i,1) - samplepoints(i);
    newband(i,2) = band(i,2) - samplepoints(i);
end
        
            
axes(axes1)
if  plotallFlag, plot(warpnocnew,tg,'MarkerSize',5); hold on; end
if ~isempty(warptest), plot(warptest(1:s_sBn) - repmat([1:s_sBn]',1,1),'kx-','MarkerSize',5); hold on;end

plot(1:s(1),newband,'r-');
plot(1:s(1),ones(s(1),1),'k-');
v = axis;
v2=max(newband(:,1));
v3=min(newband(:,2));
axis([1,s(1),min(v(3),v3)*1.15,max(v(4),v3)*1.15]);
xlabel('Sampling time (Reference)','FontSize', 12,'FontWeight','Bold');
ylabel('Sampling time (Test)','FontSize', 12,'FontWeight','Bold');
hold off;
    