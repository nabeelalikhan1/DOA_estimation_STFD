function [IF,out, peaks] = component_linking_new(tfd,orient, thr, L,theta)
global phi;
out = [];
peaks = [];
phi=theta;
%% Inputs Checkup
if(nargin == 0), fprintf(2,'ERROR (No Inputs)\n'); return;
elseif(nargin == 2)
    if(isempty(tfd)), fprintf(2,'ERROR (Input TFD is empty)\n'); return; end
    thr = 0.01; L = 32;
elseif(nargin == 3)
    if(isempty(tfd)), fprintf(2,'ERROR (Input TFD is empty)\n'); return; end
    if(isempty(thr)), thr = 0.01;
    elseif(length(thr)>1), fprintf(2,'ERROR (Threshold must be a 1x1 number)\n'); return;
    elseif(~isnumeric(thr)), fprintf(2,'ERROR (Threshold must be a 1x1 number)\n'); return;
    end
    L = 32;
elseif(nargin == 4)
    if(isempty(tfd)), fprintf(2,'ERROR (Input TFD is empty)\n'); return; end
    if(isempty(thr)), thr = 0.01;
    elseif(length(thr)>1), fprintf(2,'ERROR (Threshold must be a 1x1 number)\n'); return;
    elseif(~isnumeric(thr)), fprintf(2,'ERROR (Threshold must be a 1x1 number)\n'); return;
    end
    if(isempty(L)), L = 32;
    elseif(length(L)>1), fprintf(2,'ERROR (Threshold must be a 1x1 number)\n'); return;
    elseif(~isnumeric(L)), fprintf(2,'ERROR (Threshold must be a 1x1 number)\n'); return;
    end
end
%% main
% Create the Binary Image of Significant Local Peaks

peaks = peak_tfd(tfd, thr);
% Component Linking for IF Estimation
[IF, out] = edgelink(peaks, orient,L);
out(out > 0) = 1;
end

%% Auxiliary Functions
function peaks = peak_tfd(tfd,thresh)
% Usage:    peaks = Peak_TFD(tfd,[thresh])
%
%   Inputs: tfd ~ The time-frequency representation where local peaks are
%                 to be found
%           thresh ~ percentage of tfd maximum local peaks must be above
%                    to be considered a valid peak. Default = 0.5
%
%   Outputs: peaks ~ a binary image indicating the located peaks
%
% Notes: This function finds the local peaks in the time frequency
%        representation which are a above a predetermined threshold
%
% Written by: PhD studint/Post doc: Luke Rankine ~ January 2006
% Supervisor: Prof. B. Boashash
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    thresh = 0.5;
end
% initialize peaks
sc1 = size(tfd,1);
sc2 = size(tfd,2);
peaks = zeros(sc1,sc2);
mag_thresh = max(tfd(:))*thresh;
for i = 1:sc2
    time = tfd(:,i);
    d_time = diff(time);
    vect1 = [];
    for j = 2:length(d_time)
        if sign(d_time(j)) < 0 && sign(d_time(j-1)) > 0
            vect1 = [vect1,j];
        end
    end
    % searching for erroneous multipeaks
    if numel(vect1) ~= 0
        vect2 = vect1(1);
        close_neigh = 1; % close neighbour distance, can be varied to suit signal length
        for k = 2:length(vect1)
            if vect1(k) >= vect1(k-1) + close_neigh
                vect2 = [vect2,vect1(k)]; % otherwise forget the second point
            end
        end
    else
        vect2 = vect1;
    end
    
    for l = 1:length(vect2)
        if tfd(vect2(l),i) > mag_thresh
            peaks(vect2(l),i) = 1;
        end
    end
end
end
% EDGELINK - Link edge points in an image into lists
%
% Usage: [edgelist edgeim] = edgelink(im, minlength, location)
%
% Arguments:  im         - Binary edge image, it is assumed that edges
%                          have been thinned.
%             minlength  - Minimum edge length of interest
%             location   - Optional complex valued image holding subpixel
%                          locations of edge points. For any pixel the
%                          real part holds the subpixel row coordinate of
%                          that edge point and the imaginary part holds
%                          the column coordinate.  See NONMAXSUP.  If
%                          this argument is supplied the edgelists will
%                          be formed from the subpixel coordinates,
%                          otherwise the the integer pixel coordinates of
%                          points in 'im' are used.
%
% Returns:  edgelist - a cell array of edge lists in row,column coords in
%                      the form
%                     { [r1 c1   [r1 c1   etc }
%                        r2 c2    ...
%                        ...
%                        rN cN]   ....]
%
%           edgeim   - Image with pixels labeled with edge number.
%
%
% This function links edge points together into chains.  Where an edge
% diverges at a junction the function simply tracks one of the branches.
% The other branch is eventually processed as another edge.  These `broken
% branches' can be remerged by MERGESEG after a call to LINESEG.
%
% See also:  DRAWEDGELIST, LINESEG, MAXLINEDEV, MERGESEG, DRAWSEG, NONMAXSUP

% Acknowledgement:
% This code is inspired by David Lowe's Link.c function from the Vista image
% processing library developed at the University of British Columbia
%    http://www.cs.ubc.ca/nest/lci/vista/vista.html

% Copyright (c) 2001-2005 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% February  2001 - Original version
% September 2004 - Revised to allow subpixel edge data to be used
% January   2006 - Edited for TF Peaks by Luke Rankine

function [edgelist, edgeim] = edgelink(im, orient, minlength, location)
global EDGEIM;      % Some global variables to avoid passing (and
% copying) of arguments, this improves speed.
global ROWS;
global COLS;
EDGEIM = im ~= 0;                     % make sure image is binary.
EDGEIM = double(EDGEIM);
[ROWS, COLS] = size(EDGEIM);
edgeNo = 1;
% Perform raster scan through image looking for edge points.  When a
% point is found trackedge is called to find the rest of the edge
% points.  As it finds the points the edge image pixels are labeled
% with the -ve of their edge No
for r = 1:ROWS
    for c = 1:COLS
        if EDGEIM(r,c) == 1
            [thereIsAPoint, ~, ~] = nextpoint(r,c, orient,[]); % Find next connected
            if thereIsAPoint==0
                EDGEIM(r,c) = 0;
            end
            
        end
    end
end

for r = 1:ROWS
    for c = 1:COLS
        if EDGEIM(r,c) == 1
            edgepoints = trackedge(r,c,orient, edgeNo, minlength);
            if ~isempty(edgepoints)
                edgelist{edgeNo} = edgepoints;
                edgeNo = edgeNo + 1;
            end
        end
    end
end
if edgeNo==1
    edgelist={};
end

edgeim = -EDGEIM;  % Finally negate image to make edge encodings +ve.
% If subpixel edge locations are supplied upgrade the integer precision
% edgelists that were constructed with data from 'location'.
if nargin == 4
    for I = 1:length(edgelist)
        ind = sub2ind(size(im),edgelist{I}(:,1),edgelist{I}(:,2));
        edgelist{I}(:,1) = real(location(ind))';
        edgelist{I}(:,2) = imag(location(ind))';
    end
end

%----------------------------------------------------------------------
% TRACKEDGE
%
% Function to track all the edge points associated with a start point.  From
% a given starting point it tracks in one direction, storing the coords of
% the edge points in an array and labelling the pixels in the edge image
% with the -ve of their edge number.  When no more connected points are
% found the function returns to the start point and tracks in the opposite
% direction.  Finally a check for the overall number of edge points found is
% made and the edge is ignored if it is too short.
%
% Note that when a junction is encountered along an edge the function
% simply tracks one of the branches.  The other branch is eventually
% processed as another edge.  These `broken branches' can be remerged by
% MERGESEG.
%
% Usage:   edgepoints = trackedge(rstart, cstart, edgeNo, minlength)
%
% Arguments:   rstart, cstart   - row and column No of starting point
%              edgeNo           - the current edge number
%              minlength        - minimum length of edge to accept
%
% Returns:     edgepoints       - Nx2 array of row and col values for
%                                 each edge point.
%                                 An empty array is returned if the edge
%                                 is less than minlength long.

end
function edgepoints = trackedge(rstart, cstart, orient, edgeNo, minlength)

global EDGEIM;
edgepoints = [rstart cstart];      % Start a new list for this edge.
EDGEIM(rstart,cstart) = -edgeNo;   % Edge points in the image are
% encoded by -ve of their edgeNo.
[thereIsAPoint, r, c] = nextpoint(rstart,cstart, orient,edgepoints); % Find next connected
% edge point.
while thereIsAPoint
    
            
            edgepoints = [edgepoints             % Add point to point list
                r    c];
       
    EDGEIM(r,c) = -edgeNo;               % Update edge image
    [thereIsAPoint, r, c] = nextpoint(r,c, orient,edgepoints);
   
end
edgepoints = flipud(edgepoints);  % Reverse order of points
% Now track from original point in the opposite direction

    [thereIsAPoint, r, c] = nextpoint(rstart,cstart, orient,edgepoints);
    

while thereIsAPoint
  %  if length(edgepoints)>1
  %      if isempty(find(edgepoints(:,2)==c, 1))
            
            edgepoints = [edgepoints             % Add point to point list
                r    c];
    %    else
    %%        break;
    %    end
    %else
   %     edgepoints = [edgepoints             % Add point to point list
   %         r    c];
   % end
    EDGEIM(r,c) = -edgeNo;               % Update edge image
    
    [thereIsAPoint, r, c] = nextpoint(r,c, orient,edgepoints);
    %    if ~isempty(find(edgepoints(2:2:end)==c, 1))
    %   break;
end
% Reject short edges
Npoints = size(edgepoints,1);
if Npoints < minlength
    for i = 1:Npoints   % Clear pixels in the edge image
        EDGEIM(edgepoints(i,1), edgepoints(i,2)) = 0;
    end
    edgepoints = [];    % Return empty array
end
end

%----------------------------------------------------------------------
%
% NEXTPOINT
%
% Function finds a point that is 8 connected to an existing edge point

function [thereIsAPoint, nextr, nextc] = nextpoint(rp,cp, orient,edgepoints)
global EDGEIM;
global ROWS;
global COLS;
global phi;


% Modified by Luke Rankine, Jan. 2006
% row and column offsets for the 20(origally eight) neighbours of a
% point starting with those that are semi 4-connected.
roff = [1  0 -1  0  1  1 -1 -1  2 -2 -2  2 1  1 -1 -1 0  0 2 -2];%   3 -3 -3  3 2  2 -2 -2 1  1 3 -3  0  0 3 -3 ];
coff = [0  1  0 -1  1 -1 -1  1  1  1 -1 -1 2 -2 -2  2 2 -2 0  0]; %  2  2 -2 -2 3 -3 -3  3 3 -3 1  1  3 -3 0  0 ];
roff = [1  0 -1  0  1  1 -1 -1  2 -2 -2  2 1  1 -1 -1 0  0 2 -2  3 -3 -3  3 2  2 -2 -2 1  1 3 -3  0  0 3 -3 ];
coff = [0  1  0 -1  1 -1 -1  1  1  1 -1 -1 2 -2 -2  2 2 -2 0  0  2  2 -2 -2 3 -3 -3  3 3 -3 1  1  3 -3 0  0 ];
% 
% for i=-MM:MM
%     for j=-LL:LL
%         
%         if j==0
%             continue;
%         else
%                 rr(1,k)=i;
%                 rr(2,k)=j;
%                 rr(3,k)=i^2+j^2;
%                 k=k+1;
%         end
%         
%     end
% end

r = rp+roff;
c = cp+coff;

% Search through neighbours and return first connected edge point

for i = 1:length(r)
    if r(i) >= 1 && r(i) <= ROWS && c(i) >= 1 && c(i) <= COLS
        if EDGEIM(r(i),c(i)) == 1
            % abs(orient(r(i),c(i))-orient(rp,cp))
            diff_ang=abs(orient(r(i),c(i))-orient(rp,cp));
            if diff_ang>90
                diff_ang=180-diff_ang;
            end
            if diff_ang<phi
                if isempty(edgepoints)
                    
                    nextr = r(i);
                    nextc = c(i);
                    thereIsAPoint = 1;
                    return;
                elseif isempty(find(edgepoints(:,2)==c(i), 1))
                    
                    nextr = r(i);
                    nextc = c(i);
                    thereIsAPoint = 1;
                    return;
                end
            end
            
        end% break out and return with the data
        
    end
end
% If we get here there was no connecting next point.
nextr = 0;
nextc = 0;
thereIsAPoint = 0;
end


