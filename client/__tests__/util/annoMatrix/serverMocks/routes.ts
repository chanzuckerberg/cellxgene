import { schema } from "./schema";
import { Dataframe, KeyIndex } from "../../../../src/util/dataframe";
import { encodeMatrixFBS } from "../../../../src/util/stateManager/matrix";

const indexedSchema = {
  obsByName: Object.fromEntries(
    schema.schema.annotations.obs.columns.map((v) => [v.name, v]) ?? []
  ),
  varByName: Object.fromEntries(
    schema.schema.annotations.var.columns.map((v) => [v.name, v]) ?? []
  ),
  embByName: Object.fromEntries(
    schema.schema.layout.obs.map((v) => [v.name, v]) ?? []
  ),
};

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function makeMockColumn(s: any, length: any) {
  const { type } = s;
  switch (type) {
    case "int32":
      return new Int32Array(length).fill(Math.floor(99 * Math.random()));

    case "string":
      return new Array(length).fill("test");

    case "float32":
      return new Float32Array(length).fill(99 * Math.random());

    case "boolean":
      return new Array(length).fill(false);

    case "categorical":
      return new Array(length).fill(s.categories[0]);

    default:
      throw new Error("unkonwn type");
  }
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function getEncodedDataframe(colNames: any, length: any, colSchemas: any) {
  const colIndex = new KeyIndex(colNames);
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const columns = colSchemas.map((s: any) => makeMockColumn(s, length));
  // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
  const df = new Dataframe([length, colNames.length], columns, null, colIndex);
  const body = encodeMatrixFBS(df);
  return body;
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function dataframeResponse(colNames: any, columns: any) {
  const colIndex = new KeyIndex(colNames);
  const df = new Dataframe(
    [columns[0].length, colNames.length],
    columns,
    null,
    // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
    colIndex
  );
  const body = encodeMatrixFBS(df);
  const headers = new Headers({
    "Content-Type": "application/octet-stream",
  });
  return () => Promise.resolve({ body, init: { status: 200, headers } });
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function annotationObsResponse(request: any) {
  const url = new URL(request.url);
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const params = Array.from((url.searchParams as any).entries());
  const names = params
    // @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call.
    .filter(([k]) => k === "annotation-name")
    // @ts-expect-error ts-migrate(2345) FIXME: Argument of type '([, v]: [any, any]) => any' is n... Remove this comment to see the full error message
    .map(([, v]) => v);
  // @ts-expect-error ts-migrate(2538) FIXME: Type 'unknown' cannot be used as an index type.
  if (!names.every((n) => indexedSchema.obsByName[n])) {
    return Promise.reject(new Error("bad obs annotation name in URL"));
  }
  // @ts-expect-error ts-migrate(2538) FIXME: Type 'unknown' cannot be used as an index type.
  const colSchemas = names.map((n) => indexedSchema.obsByName[n]);
  const body = getEncodedDataframe(
    names,
    schema.schema.dataframe.nObs,
    colSchemas
  );

  const headers = new Headers({
    "Content-Type": "application/octet-stream",
  });
  return Promise.resolve({
    body,
    init: { status: 200, headers },
  });
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function annotationVarResponse(request: any) {
  const url = new URL(request.url);
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const params = Array.from((url.searchParams as any).entries());
  const names = params
    // @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call.
    .filter(([k]) => k === "annotation-name")
    // @ts-expect-error ts-migrate(2345) FIXME: Argument of type '([, v]: [any, any]) => any' is n... Remove this comment to see the full error message
    .map(([, v]) => v);
  // @ts-expect-error ts-migrate(2538) FIXME: Type 'unknown' cannot be used as an index type.
  if (!names.every((n) => indexedSchema.varByName[n])) {
    return Promise.reject(new Error("bad var annotation name in URL"));
  }
  // @ts-expect-error ts-migrate(2538) FIXME: Type 'unknown' cannot be used as an index type.
  const colSchemas = names.map((n) => indexedSchema.varByName[n]);
  const body = getEncodedDataframe(
    names,
    schema.schema.dataframe.nVar,
    colSchemas
  );

  const headers = new Headers({
    "Content-Type": "application/octet-stream",
  });
  return Promise.resolve({
    body,
    init: { status: 200, headers },
  });
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function layoutObsResponse(request: any) {
  const url = new URL(request.url);
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const params = Array.from((url.searchParams as any).entries());
  // @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call.
  const names = params.filter(([k]) => k === "layout-name").map(([, v]) => v);
  // @ts-expect-error ts-migrate(2538) FIXME: Type 'unknown' cannot be used as an index type.
  if (!names.every((n) => indexedSchema.embByName[n])) {
    return Promise.reject(new Error("bad layout name in URL"));
  }
  // @ts-expect-error ts-migrate(2538) FIXME: Type 'unknown' cannot be used as an index type.
  const dims = names.map((n) => indexedSchema.embByName[n].dims).flat();
  const colSchemas = names
    // @ts-expect-error ts-migrate(2538) FIXME: Type 'unknown' cannot be used as an index type.
    .map((n) => [indexedSchema.embByName[n], indexedSchema.embByName[n]])
    .flat();
  const body = getEncodedDataframe(
    dims,
    schema.schema.dataframe.nObs,
    colSchemas
  );

  const headers = new Headers({
    "Content-Type": "application/octet-stream",
  });
  return Promise.resolve({
    body,
    init: { status: 200, headers },
  });
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function dataVarResponse(request: any) {
  const url = new URL(request.url);
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const params = Array.from((url.searchParams as any).entries());

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const colNames = params.map((v) => `${(v as any)[0]}/${(v as any)[1]}`);
  const colSchemas = colNames.map(() => schema.schema.dataframe);
  const body = getEncodedDataframe(
    colNames,
    schema.schema.dataframe.nObs,
    colSchemas
  );

  const headers = new Headers({
    "Content-Type": "application/octet-stream",
  });
  return Promise.resolve({
    body,
    init: { status: 200, headers },
  });
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function responder(request: any) {
  const url = new URL(request.url);
  const { pathname } = url;
  if (pathname.endsWith("/annotations/obs")) {
    return annotationObsResponse(request);
  }
  if (pathname.endsWith("/annotations/var")) {
    return annotationVarResponse(request);
  }
  if (pathname.endsWith("/layout/obs")) {
    return layoutObsResponse(request);
  }
  if (pathname.endsWith("/data/var")) {
    return dataVarResponse(request);
  }
  return Promise.reject(new Error("bad URL"));
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function withExpected(expectedURL: any, expectedParams: any) {
  /*
  Do some additional error checking
  */
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  return (request: any) => {
    // if URL is bogus, reject the promise
    const url = new URL(request.url);
    if (!url.pathname.endsWith(expectedURL)) {
      return Promise.reject(new Error("Unexpected URL!"));
    }
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const params = Array.from((url.searchParams as any).entries()).sort(
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type '(a: unknown, b: unknown) => bool... Remove this comment to see the full error message
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (a, b) => (a as any)[0] < (b as any)[0]
    );
    expectedParams = expectedParams
      .slice() // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      .sort((a: any, b: any) => a[0] < b[0]);

    if (
      params.length !== expectedParams.length ||
      !params.every(
        (p, i) =>
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          (p as any)[0] === expectedParams[i][0] && // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          (p as any)[1] === expectedParams[i][1]
      )
    ) {
      return Promise.reject(new Error("unexpected name requested in URL"));
    }

    return responder(request);
  };
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function annotationsObs(names: any) {
  return withExpected(
    "/annotations/obs",
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    names.map((name: any) => ["annotation-name", name])
  );
}
