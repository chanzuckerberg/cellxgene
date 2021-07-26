/*
Define globals which are present in the client, but not in node (and therefore not in
the jest test environment).
*/

import { TextDecoder, TextEncoder } from "util";

// @ts-expect-error ts-migrate(2322) FIXME: Type 'typeof TextDecoder' is not assignable to typ... Remove this comment to see the full error message
global.TextDecoder = TextDecoder;
global.TextEncoder = TextEncoder;
