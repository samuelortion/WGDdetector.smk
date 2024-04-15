#!/usr/bin/env -S awk -f

/^[^\#]/ && $10 < threshold {
    print $0
}