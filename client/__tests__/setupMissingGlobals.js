/*
Define globals which are present in the client, but not in node (and therefore not in
the jest test environment).
*/

import { TextDecoder, TextEncoder } from "util";

global.TextDecoder = TextDecoder;
global.TextEncoder = TextEncoder;
