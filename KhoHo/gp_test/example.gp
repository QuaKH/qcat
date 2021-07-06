/*
    Determine if integer is a prime
*/

is_prime(n) = {
    local(dividend);
    dividend = 2;
    while (
        dividend <= sqrt(n),
        if (
            n % dividend == 0,
            return(0),
            dividend = dividend + 1
        )
    );
    return(1);
}