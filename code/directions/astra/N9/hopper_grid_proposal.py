"""Proposed bounded matrix study for Hopper. Not submitted by this review.

Example: python3 -I hopper_grid_proposal.py --kind dh --a 3.715 --denom 1000
Use an environment containing python-flint, numpy, scipy and mpmath.
"""
import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from dh_certify import certify_matrix
from dh_arb import covariance as dh_covariance
from zeta_arb import psi_arb, covariance as zeta_covariance
from flint import ctx

parser=argparse.ArgumentParser()
parser.add_argument('--kind',choices=['zeta','dh'],required=True)
parser.add_argument('--a',required=True)
parser.add_argument('--denom',type=int,required=True)
parser.add_argument('--bits',type=int,default=256)
parser.add_argument('--output',default='hopper_grid_result.json')
args=parser.parse_args()
n=Fraction(args.a)*args.denom
if n.denominator!=1:
    parser.error('a must be a multiple of 1/denom')
n=int(n)
if not 1<=n<=8000:
    parser.error('This proposal is limited to 8000 increments')
ctx.prec=args.bits
if args.kind=='zeta':
    c=zeta_covariance(psi_arb(args.denom,n))
else:
    # The proposed denominators 400 and 1000 have terminating decimals.
    c=dh_covariance(str(1/args.denom),n,args.bits)
result={'kind':args.kind,'a':args.a,'denom':args.denom,'bits':args.bits,
        'certificate':certify_matrix(c)}
Path(args.output).write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
