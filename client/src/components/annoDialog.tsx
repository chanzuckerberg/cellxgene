import React from "react";
import { Button, Tooltip, Dialog, Classes, Colors } from "@blueprintjs/core";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
type State = any;

// eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
class AnnoDialog extends React.PureComponent<{}, State> {
  // eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
  constructor(props: {}) {
    super(props);
    this.state = {};
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isActive' does not exist on type 'Readon... Remove this comment to see the full error message
      isActive,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'text' does not exist on type 'Readonly<{... Remove this comment to see the full error message
      text,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'title' does not exist on type 'Readonly<... Remove this comment to see the full error message
      title,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'instruction' does not exist on type 'Rea... Remove this comment to see the full error message
      instruction,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'cancelTooltipContent' does not exist on ... Remove this comment to see the full error message
      cancelTooltipContent,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'errorMessage' does not exist on type 'Re... Remove this comment to see the full error message
      errorMessage,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'validationError' does not exist on type ... Remove this comment to see the full error message
      validationError,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoSelect' does not exist on type 'Read... Remove this comment to see the full error message
      annoSelect,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoInput' does not exist on type 'Reado... Remove this comment to see the full error message
      annoInput,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'secondaryInstructions' does not exist on... Remove this comment to see the full error message
      secondaryInstructions,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'secondaryInput' does not exist on type '... Remove this comment to see the full error message
      secondaryInput,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'handleCancel' does not exist on type 'Re... Remove this comment to see the full error message
      handleCancel,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'handleSubmit' does not exist on type 'Re... Remove this comment to see the full error message
      handleSubmit,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'primaryButtonText' does not exist on typ... Remove this comment to see the full error message
      primaryButtonText,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'secondaryButtonText' does not exist on t... Remove this comment to see the full error message
      secondaryButtonText,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'handleSecondaryButtonSubmit' does not ex... Remove this comment to see the full error message
      handleSecondaryButtonSubmit,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'primaryButtonProps' does not exist on ty... Remove this comment to see the full error message
      primaryButtonProps,
    } = this.props;

    return (
      <Dialog icon="tag" title={title} isOpen={isActive} onClose={handleCancel}>
        <form
          onSubmit={(e) => {
            e.preventDefault();
          }}
        >
          <div className={Classes.DIALOG_BODY}>
            <div style={{ marginBottom: 20 }}>
              <p>{instruction}</p>
              {annoInput || null}
              <p
                style={{
                  marginTop: 7,
                  visibility: validationError ? "visible" : "hidden",
                  color: Colors.ORANGE3,
                }}
              >
                {errorMessage}
              </p>
              {/* we might rename, secondary button and secondary input are not related */}
              {secondaryInstructions && (
                <p style={{ marginTop: secondaryInstructions ? 20 : 0 }}>
                  {secondaryInstructions}
                </p>
              )}
              {secondaryInput || null}
            </div>
            {annoSelect || null}
          </div>
          <div className={Classes.DIALOG_FOOTER}>
            <div className={Classes.DIALOG_FOOTER_ACTIONS}>
              <Tooltip content={cancelTooltipContent}>
                <Button onClick={handleCancel}>Cancel</Button>
              </Tooltip>
              {/* we might rename, secondary button and secondary input are not related */}
              {handleSecondaryButtonSubmit && secondaryButtonText ? (
                <Button
                  onClick={handleSecondaryButtonSubmit}
                  disabled={!text || validationError}
                  intent="none"
                  type="button"
                >
                  {secondaryButtonText}
                </Button>
              ) : null}
              <Button
                {...primaryButtonProps} // eslint-disable-line react/jsx-props-no-spreading -- Spreading props allows for modularity
                onClick={handleSubmit}
                disabled={!text || validationError}
                intent="primary"
                type="submit"
              >
                {primaryButtonText}
              </Button>
            </div>
          </div>
        </form>
      </Dialog>
    );
  }
}

export default AnnoDialog;
