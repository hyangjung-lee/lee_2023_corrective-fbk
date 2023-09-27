
function [] = save_figureToPDF( handle, file_baseName, fig_area, paper_DPI )

% Usage: save_figureToPDF( handle, fileBaseName, figArea, paper_DPI)
%
% Saves a given figure specified by 'handle' into PDF format file with the name 'fileBaseName'.
% You can also specify the area of figure which is converted to a PDF file.
% If you want to adjust the resolution, it is possible to set the pixes per inch in PDF file.
%
% handle:        Handle to a Matlab figure
%
% fileBaseName:  The name of PDF file, not including '.pdf' extension.
%
% 'figArea':     (optional) The area to be saved into PDF file.
%                Its format is [bottom left width height].
%                The unit is a 'normalized' unit which ranges from 0 to 1.
%                If this parameter isn't specified,
%                the whole area of figure will be printed to PDF file.
%
% 'paper_DPI':   (optional) Pixels per inch in PDF file.
%                If this value isn't specified, DPI of screen will be used instead.
%
% EXAMPLE:
%    save_figureToPDF( handle, 'figures/plot1' );
%    save_figureToPDF( handle, 'fig_1', [.2 .1 .8 .9] );
%    save_figureToPDF( handle, 'fig_2', [.2 .1 .8 .9], 200 );
%
% Copyright (C) Daeseob Lim (daeseob at gmail dot com), 2013.

	if ( ~exist( 'fig_area', 'var' ) || isempty(fig_area) )
		fig_area = [0 0 1 1];
	end

	pixel_pos = getpixelposition( handle );
	% NOTE) If we cut the figure according to 'fig_area',
	%	the right and the bottom portions of the figure are fully shown.
	%	So, we give a little margin to right and bottom areas.
	deltaX_perPixel	= 1 / pixel_pos(3);
	deltaY_perPixel	= 1 / pixel_pos(4);
	right_margin	= 2;	% in pixels
	bottom_margin	= 2;	% in pixels

	h_pos	= fig_area(1);
	v_pos	= fig_area(2) - bottom_margin * deltaY_perPixel;
	width	= fig_area(3) + right_margin * deltaX_perPixel;
	height	= fig_area(4) + bottom_margin * deltaY_perPixel;

	% Determine the size of output image in terms of paper size (in cm)
	screen_DPI = get( 0, 'ScreenPixelsPerInch' );
	wholeFigureSize	= pixel_pos(3:4) / screen_DPI * 2.54;					% in cm
	validFigureSize	= wholeFigureSize .* [width height];	% in cm

	paperSize		= validFigureSize;
	leftMargin		= wholeFigureSize(1) * h_pos;
	bottomMargin	= wholeFigureSize(2) * v_pos;
	paperPosition	= [-leftMargin, -bottomMargin, wholeFigureSize];

	% Set the proper properties
	set( handle, 'PaperUnits', 'centimeters' );
	set( handle, 'PaperSize', paperSize );
	set( handle, 'PaperPosition', paperPosition );
	set( handle, 'InvertHardcopy', 'off' );
	
	% Print out to a PDF file (NOTE: Without '-r' option, default DPI is set to 600.)
	fileFullName = [file_baseName '.pdf'];
	if ( exist('paper_DPI','var') )
		print( handle, fileFullName, '-dpdf', ['-r' num2str(paper_DPI)] );
	else
		default_DPI = 600;
		print( handle, fileFullName, '-dpdf', ['-r' num2str(default_DPI)] );
	end

end
