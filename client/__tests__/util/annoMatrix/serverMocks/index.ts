export const baseDataURL = "https://a.fake.url/api/v0.2";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
(window as any).CELLXGENE = {
  API: {
    prefix: baseDataURL,
    version: "v0.2/",
  },
};

export { schema } from "./schema";
export * from "./routes";
