function str = random_string(length, type)
    %RANDOM_STRING - Generate a random string
    %
    %When a name is needed but none is provided, never fear
    %RANDOM_STRING is here!  As configured returns a random
    %set of upper case letters (default of 10 long)
    %
    % Input:
    % @length
    %   value   - Length of the random string to generate
    %   default - 10
    %   type    - int
    % @type
    %   value   - Type of characters to include.
    %             1 = A-Z, 2 = a-z, 3=letters/#'s/special chars
    %   default - 1
    %   type    - int
    %
    % Return:
    % @str
    %   value - The generated random string
    %   type  - str

    % If no length is specified, use 10
    if (nargin < 1)
        length = 10;
    else
        if (length < 1)
            error(sprintf('"length" argument must be positive in %s',mfilename()))
        elseif (mod(length,1) ~= 0)
            error(sprintf('"length" argument must be an integer in %s',mfilename()))
        end
    end
    
    % If the type has not been specified,
    % use upper case letters
    if (nargin < 2) type = 1; end
    
    %ASCII Character codes
    if (type == 1)
        % Upper case letters
        % A = 65 ... Z = 90
        str = char(65 + floor(26 .* rand(length,1)))';
    elseif (type == 2)
        % Lower case letters
        % a = 97 ... z = 122
        str = char(97 + floor(26 .* rand(length,1)))';
    elseif (type == 3)
        % Upper case, lower case, numbers, and special characters
        % !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
        str = char(33 + floor(94 .* rand(length,1)))';
    else
        error(sprintf('Unrecognized "type" argument in %s',mfilename()))
    end
end
