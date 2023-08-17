function AB = resex_mc_correlate(A,B,dt)
%calculates and returns discrete correlation between A and B
%assumes A and B are row/column vectors
% tic
if ~isrow(A)||~isrow(B); error("Input has wrong dimensions"); end %Important! The faster xcorr method below only works if A and B are rows
n = numel(A);
% m = numel(B);
% if m~=n; error('Can not handle this'); end
% if size(B,1) ~= size(A,1); B = B'; end
% % 
% AB = zeros(size(A));
% for c_n = 0:n-1
%     for c_m = 1:n-c_n
%        AB(c_n+1) = AB(c_n+1) + (A(c_m)*B(c_m+c_n))*dt; 
%     end
% end

% % 
AB  = xcorr(A,B)*dt; 
AB = AB(1:n);
AB = fliplr(AB);
% % toc
end