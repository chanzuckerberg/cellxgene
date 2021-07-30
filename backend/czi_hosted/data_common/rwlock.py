# -*- coding: utf-8 -*-
""" rwlock.py

    A class to implement read-write locks on top of the standard threading
    library.

    This is implemented with two mutexes (threading.Lock instances) as per this
    wikipedia pseudocode:

    https://en.wikipedia.org/wiki/Readers%E2%80%93writer_lock#Using_two_mutexes

    Code written by Tyler Neylon at Unbox Research.

    This file is public domain.

    Modified to add a w_demote function to convert a writer lock to a reader lock
"""


# _______________________________________________________________________
# Imports

from contextlib import contextmanager
from threading import Lock


# _______________________________________________________________________
# Class


class RWLock(object):
    """RWLock class; this is meant to allow an object to be read from by
    multiple threads, but only written to by a single thread at a time. See:
    https://en.wikipedia.org/wiki/Readers%E2%80%93writer_lock

    Usage:

        from rwlock import RWLock

        my_obj_rwlock = RWLock()

        # When reading from my_obj:
        with my_obj_rwlock.r_locked():
            do_read_only_things_with(my_obj)

        # When writing to my_obj:
        with my_obj_rwlock.w_locked():
            mutate(my_obj)
    """

    def __init__(self):

        self.w_lock = Lock()
        self.num_r_lock = Lock()
        self.num_r = 0

        # The d_lock is needed to handle the demotion case,
        # so that the writer can become a reader without releasing the w_lock.
        # the d_lock is held by the writer, and prevents any other thread from taking the
        # num_r_lock during that time, which means the writer thread is able to take the
        # num_r_lock to update the num_r.
        self.d_lock = Lock()

    # ___________________________________________________________________
    # Reading methods.

    def r_acquire(self):
        self.d_lock.acquire()
        self.num_r_lock.acquire()
        self.num_r += 1

        if self.num_r == 1:
            self.w_lock.acquire()

        self.num_r_lock.release()
        self.d_lock.release()

    def r_release(self):
        assert self.num_r > 0
        self.num_r_lock.acquire()
        self.num_r -= 1
        if self.num_r == 0:
            self.w_lock.release()

        self.num_r_lock.release()

    @contextmanager
    def r_locked(self):
        """This method is designed to be used via the `with` statement."""
        try:
            self.r_acquire()
            yield
        finally:
            self.r_release()

    # ___________________________________________________________________
    # Writing methods.

    def w_acquire(self):
        self.d_lock.acquire()
        self.w_lock.acquire()

    def w_acquire_non_blocking(self):
        # if d_lock and w_lock can be acquired without blocking, acquire and return True,
        # else immediately return False.
        if self.d_lock.acquire(blocking=False):
            if self.w_lock.acquire(blocking=False):
                return True
            else:
                self.d_lock.release()
        return False

    def w_release(self):
        self.w_lock.release()
        self.d_lock.release()

    def w_demote(self):
        """demote a writer lock to a reader lock"""

        # the d_lock is already held from w_acquire.
        # releasing the d_lock at the end of this function allows multiple readers.
        # incrementing num_r makes this thread one of those readers.
        self.num_r_lock.acquire()
        self.num_r += 1
        self.num_r_lock.release()
        self.d_lock.release()

    @contextmanager
    def w_locked(self):
        """This method is designed to be used via the `with` statement."""
        try:
            self.w_acquire()
            yield
        finally:
            self.w_release()
