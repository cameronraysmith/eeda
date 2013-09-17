function s = trapezoid(t,m)
%TREPEZOID Trapezoid wave generation.
%   TRAPEZOID(T) generates a square wave with period 2*Pi for the
%   elements of time vector T.  TRAPEZOID(T) is like SQUARE(T), only
%   it creates a square wave with peaks of +1 to -1 instead of
%   a sine wave.
%   
%   For example, generate a trapezoidal wave with 24 time unit period
%   and rise/fall time equal to 24*pi/8/(2*pi):
%        t=linspace(0,72,10000);
%        y=trapezoid(t*2*pi/24,8); plot(t,y,'r.');
%
%   See also SQUARE

%   Author(s): Jordan Firth, 2 Sep 2010
%   $Revision: 0.0 $  $Date: 2007/12/14 15:06:21 $

tmp = mod(t+pi/m,2*pi);

a = (tmp < pi/m);
b = (tmp >= pi/m & tmp < pi);
c = (tmp >= pi & tmp < pi+pi/m);

rise = m*tmp/pi;
fall = -m*(tmp-pi)/pi+1;
nodd = a.*rise + b + c.*fall;

s = 2*nodd-1;