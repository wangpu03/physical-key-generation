in this document, we simulate the OFDM communication, ambient backscatter communication between two BDs.
1. ofdm_sim
    the OFDM transmitting and receiving with different modulation
2. ambc_ofdm
    the direct transmitting with one channel over OFDM
    the undirect transmitting with two combination channel over OFDM
    discuss the multiplication of these two direct and undirect informaiton 
3. ambc_ofdm_02
    the completely ambient backscatter communication over OFDM of the paper system model 
    direct signal, backscatter signal, combination signals
    calculate the direct and backscatter signal from the combination signals. 
    with same transmission information, the power of multiplication of two signals only consisting of CP part are the same.
4. ambc_ofdm_03
    based on the ambc_ofdm_02.m, we add some noise and other parameters, such as the backscatter coefficient, different path gain in different channel.