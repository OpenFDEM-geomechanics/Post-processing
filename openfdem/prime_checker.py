import concurrent.futures
import math 

def is_prime(n):
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False

    sqrt_n = int(math.floor(math.sqrt(n)))
    for i in range(3, sqrt_n + 1, 2):
        if n % i == 0:
            return False
    return True

def run_checker(PRIMES):
    results =[]

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = executor.map(is_prime, PRIMES)
        for number, prime in zip(PRIMES, futures):
            print('%d is prime: %s' % (number, prime))
            results.append(number)

    print(results)