import { Position, Toaster, Intent } from "@blueprintjs/core";

/* styles */
// @ts-expect-error ts-migrate(2307) FIXME: Cannot find module './menubar.css' or its correspo... Remove this comment to see the full error message
import styles from "./toasters.css";

/** Singleton toaster instance. Create separate instances for different options. */

const ToastTopCenter = Toaster.create({
  className: "recipe-toaster",
  position: Position.TOP,
  maxToasts: 4,
});

/*
A "user" error - eg, bad input
*/
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export const postUserErrorToast = (message: any) =>
  ToastTopCenter.show({ message, intent: Intent.WARNING });

/*
A toast the user must dismiss manually, because they need to act on its information,
ie., 8 bulk add genes out of 40 were bad. Manually see which ones and fix.
*/
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export const keepAroundErrorToast = (message: any) =>
  ToastTopCenter.show({ message, timeout: 0, intent: Intent.WARNING });

/*
a hard network error
*/
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export const postNetworkErrorToast = (message: any, key = undefined) =>
  ToastTopCenter.show(
    {
      message,
      timeout: 30000,
      intent: Intent.DANGER,
    },
    key
  );

/*
Async message to user
*/
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export const postAsyncSuccessToast = (message: any) =>
  ToastTopCenter.show({
    message,
    timeout: 10000,
    intent: Intent.SUCCESS,
  });

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export const postAsyncFailureToast = (message: any) =>
  ToastTopCenter.show({
    message,
    timeout: 10000,
    intent: Intent.WARNING,
  });

/*
Dataset opened in new tab
 */
export const postExplainNewTab = (message: string) => {
  ToastTopCenter.show({
    className: styles.newTabToast,
    message,
    timeout: 5000,
    intent: Intent.PRIMARY,
  });
};
