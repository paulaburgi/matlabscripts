function oops(n)
%
% OOPS Delete the last object plotted on the axes.
%      Repeating "oops" erases further back in time.
%      OOPS does not work for title and labels; 
%      to erase these, use "title('')" or "xlabel('')"
%
%      From Knight, A., "Basics of Matlab and Beyond", p. 113
if nargin == 0
	n = 1;
end
h = get(gca,'children');
delete(h(1:n));