function str = TxRx(num)
    b = floor((num-1)/8); 
    if b == 0
        band = '900 MHz';
    else 
        band = '2 GHz';
    end 
    num = num - b * 8; 
    n = ceil(num/2);

    if mod(num, 2) == 1
       type = 'Amp';
    else
        type = 'Phs';
    end
    
    tx = mod(n-1,2)+1;
    rx = ceil(n/2);
    str = strcat('Tx', num2str(tx), 'Rx', num2str(rx), type, '_', band);
end